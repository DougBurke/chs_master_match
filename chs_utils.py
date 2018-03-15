"""
Utility routines.
"""

import glob
import os

import numpy as np

import pycrates
import region

from coords.format import deg2ra, deg2dec


def find_single_match(pat):
    """Return the file name matching the pattern.

    If there's a single match, return it, otherwise file (so
    if there's no matches or multiple matches).

    Parameters
    ----------
    pat : str
        The pattern passed to glob.glob.

    Returns
    -------
    filename : str
        The matching file name.

    """

    ms = glob.glob(pat)
    if len(ms) == 1:
        return ms[0]
    elif len(ms) == 0:
        raise IOError("No match for [{}]".format(pat))
    else:
        raise IOError("Multiple matches to [{}]\n{}".format(pat,
                                                            ms))


def find_mrgsrc3(stack, mrgsrc3dir):
    """Return the mrgsrc3 file name.

    Parameters
    ----------
    stack : str
        The stack name.
    mrgsrc3dir : str
        The location of the mrgsrc3 files.

    """

    pat = "{}*mrgsrc3.fits*".format(stack)
    pat = os.path.join(mrgsrc3dir, pat)
    return find_single_match(pat)


def find_stkfov3(stack, fov3dir):
    """Return the fov3 file name.

    Parameters
    ----------
    stack : str
        The stack name.
    fov3dir : str
        The location of the fov3 files.

    """

    pat = "{}*fov3.fits*".format(stack)
    pat = os.path.join(fov3dir, pat)
    return find_single_match(pat)


def find_stkevt3(stack, evt3dir):
    """Return the evt3 file name.

    Parameters
    ----------
    stack : str
        The stack name.
    evt3dir : str
        The location of the evt3 files.

    """

    pat = "{}*evt3.fits*".format(stack)
    pat = os.path.join(evt3dir, pat)
    return find_single_match(pat)


# I do have code to do this in a "nicer" way, in that it uses the
# C API of the region library to make constructing and passing
# around regions a lot nicer in Python. Alternatively, the
# new-to-CIAO 4.10 region API could be used. However, it should
# not be needed here (the regions being created are expected to
# be significantly larger than one pixel).
#
def make_region_string(xy, ndp=2):
    """Create a region string representing the polygon data.

    Parameters
    ----------
    xy : numpy_array
        The vertexes of the polygon. No check is made to see
        if the list is closed, assuming that the region library
        handles this, but any non-finite values are removed.
        The shape is (2, npts) and the x and y values for a point
        are assumed to either both be good (finite) or bad (not
        finite), as well as there being at least three valid points.
    ndp : int, optional
        The number of decimal places to report each position
        as; if None then let Python decide.

    Returns
    -------
    regstr : str
        The polygon region as a string.

    Notes
    -----
    To simplify this the coordinates are converted to 2 decimal places
    (by default). Other options are to just let Python calculate the
    number of decimal places, or to use a lower-level interface to the
    region library where arrays can be sent directly. For the
    convex-hull case this should not cause problems (since we expect
    the hulls to be significantly larger than a sky pixel).

    """

    if xy.ndim != 2 or xy.shape[0] != 2 or xy.shape[1] < 3:
        raise ValueError("unxepected shape: {}".format(xy.shape))

    x = xy[0]
    y = xy[1]

    xidx = np.isfinite(x)
    yidx = np.isfinite(y)
    if (xidx != yidx).any():
        raise ValueError("Non-finite issue: {}".format(xy))

    x = x[xidx]
    y = y[xidx]

    if ndp is None:
        fmt = "{},{}"
    else:
        fmt = "{{:.{0}f}},{{:.{0}f}}".format(ndp)

    coords = [fmt.format(*z) for z in zip(x, y)]
    return "polygon({})".format(",".join(coords))


def check_overlap(reg1, reg2):
    """Do the two regions overlap?

    Parameters
    ----------
    reg1, reg2 : str
        The output of make_region_string

    Returns
    -------
    overlap : bool
        True if there's an overlap.

    Notes
    -----
    Using a lower-level interface to the region library would allow
    direct use of region objects, rather than passing around strings.

    This returns true if the overlap area is not 0. This means that
    even the smallest of overlaps (down to the accuracy of the
    intersection) is included.

    For now a "fast" mode to increase the overlap-area checks is not
    used. Perhaps it should be - or GPC should be used.
    """

    rstr = "{}&{}".format(reg1, reg2)
    reg = region.regParse(rstr)

    # Could check the bounds and if large, switch to a larger binning
    # value.
    area = region.regArea(reg)
    return area > 0.0


def validate_polygon(poly, report=True, label=None):
    """Return a closed polygon.

    A screen message is displayed if the polygon is not closed,
    unless ``report`` is not set.

    Parameters
    ----------
    poly : NumPy array
        The shape must be (2, npts), with npts > 2.
    report : bool, optional
        If `True` then a message is displayed to stdout if the polygon
        fails a validation check (presently only if it is not
    label : str or None, optional
        Added to the end of the WARNING message if a problem is
        reported.

    Returns
    -------
    validated : NumPy array
        This is either (2, n) where n is normally npts or npts + 1,
        but could be less when the input contains non-finite values.

    Raises
    ------
    ValueError
        If the non-finite points do not match for the X and Y axes.

    """

    assert poly.ndim == 2
    assert poly.shape[0] == 2
    assert poly.shape[1] > 2

    x = poly[0]
    y = poly[1]

    xidx = np.isfinite(x)
    yidx = np.isfinite(y)
    if (xidx != yidx).any():
        raise ValueError("Non-finite values do not match")

    x = x[xidx]
    y = y[xidx]

    if x[0] != x[-1] or y[0] != y[-1]:
        if report:
            msg = "WARNING: input polygon is not closed"
            if label is not None:
                msg += " " + label
            print(msg)

        x = np.append(x, x[0])
        y = np.append(y, y[0])

    return np.vstack((x, y))


def make_spatial_filter_range(rmin, rmax, dmin, dmax):
    """Return a spatial filter representing this box.

    Parameters
    ----------
    rmin, rmax : number
        Minimum and maximum Right Ascension, in decimal degrees.
    dmin, dmax : number
        Minimum and maximum Declination, in decimal degrees.

    Returns
    -------
    filter_expr : str
        A Data Model spatial filter for the SKY coordinate
        system (the limits are given as sexagessmial values).

    """

    if rmin >= rmax:
        raise ValueError("rmax must be > rmin")
    if dmin >= dmax:
        raise ValueError("dmax must be > dmin")

    # filter in celestial coordinates
    fmt = ':'
    rmin_s = deg2ra(rmin, fmt)
    rmax_s = deg2ra(rmax, fmt)
    dmin_s = deg2dec(dmin, fmt)
    dmax_s = deg2dec(dmax, fmt)

    # want to keep the precision down, so drop all but the
    # second decimal point, and do not do proper rounding
    #
    def cleanval(inval):
        pos = inval.find('.')
        if pos < 0:
            return inval

        if inval[pos + 1:].find('.') != -1:
            raise ValueError("Multiple . in pos: {}".format(inval))
        return inval[:pos + 3]

    rmin_s = cleanval(rmin_s)
    rmax_s = cleanval(rmax_s)
    dmin_s = cleanval(dmin_s)
    dmax_s = cleanval(dmax_s)

    return '[sky=RECT({},{},{},{})]'.format(rmin_s, dmin_s,
                                            rmax_s, dmax_s)


def read_hulls_from_mrgsrc3(mrgsrc3):
    """Return the HULL info.

    Parameters
    ----------
    mrgsrc3 : str
        The name of the mrgsrc3 file.

    Returns
    -------
    hulls : list of dict
        The hull coordinates and metadata. It can be empty. Each hull
        is represented by a dict with keys: stack, component, band,
        mancode (True if mancode is not 0), pos, eqpos, and transform
        (which contains the WCS). The polygons are closed and only
        contain finite values.

    """

    infile = "{}[MEXTSRC][status=0]".format(mrgsrc3)
    cr = pycrates.read_file(infile)
    nr = cr.get_nrows()
    if nr == 0:
        return []

    stack = cr.get_key_value('STACK_ID')
    assert stack is not None

    tr = cr.get_transform('EQSRC').copy()

    # The manual code is converted from a bit to a simple value
    # (not in a way that encodes the flag settings). Treat as
    # a boolean (0 vs non-zero)
    #
    zs = zip(cr.COMPONENT.values,
             cr.EBAND.values,
             cr.POS.values,
             cr.get_column('EQSRC').values,
             cr.MAN_CODE.values.sum(axis=1)
             )

    out = []
    for cpt, eband, pos, eqsrc, mancode in zs:

        pos = validate_polygon(pos, report=False)
        eqsrc = validate_polygon(eqsrc, report=False)
        assert pos.shape == eqsrc.shape

        out.append({'stack': stack,
                    'component': cpt,
                    # I think there may still be some upper case band values
                    'band': eband.lower(),
                    'mancode': mancode > 0,
                    'pos': pos.copy(),
                    'eqpos': eqsrc.copy(),
                    'transform': tr,
                    })

    return out


def read_master_hulls(chsfile):
    """Read in hull information.

    Parameters
    ----------
    chsfile : str
        The name of the file created by chs_create_initial_masters.py

    Returns
    -------
    hullmatch, hulllist, metadata : dict, dict, dict
        The contents of the HULLMATCH and HULLLIST block, and metadata
        about the file (e.g. ensemble and CHSVER value). The
        hullmatch and hulllist dicts have keys of Master_Id.
        If HULLMATCH and HULLLIST can not be found then SRCMATCH
        and SRCLIST are used.
    """

    ds = pycrates.CrateDataset(chsfile, mode='r')

    prefix = 'HULL'
    try:
        cr = ds.get_crate('{}MATCH'.format(prefix))
    except IndexError:
        prefix = 'SRC'
        cr = ds.get_crate('{}MATCH'.format(prefix))

    hullmatch = {}

    # Map from stack id to the "number" of the stack in the ensemble
    ensemblemap = {}
    for i in range(cr.get_key_value('STKIDNUM')):
        key = 'STKID{:03d}'.format(i)
        stkid = cr.get_key_value(key)
        assert stkid not in ensemblemap, stkid
        ensemblemap[stkid] = i

    zs = zip(cr.Master_Id.values,
             cr.STACKID.values,
             cr.COMPONENT.values,
             cr.Match_Type.values)
    for mid, stackid, component, mtype in zs:
        store = {'master_id': mid,
                 'stack': stackid,
                 'component': component,
                 'match_type': mtype}
        try:
            hullmatch[mid].append(store)
        except KeyError:
            hullmatch[mid] = [store]

    cr = ds.get_crate('{}LIST'.format(prefix))

    hulllist = {}
    zs = zip(cr.Master_Id.values,
             cr.STATUS.values,
             cr.BASE_STK.values,
             cr.NVERTEX.values,
             cr.EQPOS.values)
    for mid, status, base_stk, nvertex, eqpos in zs:
        assert mid not in hulllist

        hulllist[mid] = {'master_id': mid,
                         'status': status,
                         'base_stk': base_stk,
                         'eqpos': eqpos[:, :nvertex]}

    metadata = {'ensemble': cr.get_key_value('ENSEMBLE'),
                'ensemblemap': ensemblemap,
                'revision': cr.get_key_value('CHSVER')}

    return hullmatch, hulllist, metadata
