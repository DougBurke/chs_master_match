"""
Utility routines.
"""

import glob
import os
import sys

import numpy as np

import pycrates
import region

from coords.format import deg2ra, deg2dec


_logmsg = set([])


def logmsg(msg):
    """Display the message to STDOUT if it hasnt' already been reported.

    The message is preceeded by 'LOG: '.

    Parameters
    ----------
    msg : str
        The message to display (actually can be anything that can
        be converted to a string and has an equality test).

    """

    if msg in _logmsg:
        return

    print("LOG: {}".format(msg))
    _logmsg.add(msg)


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


def get_revision(crate):
    """Return the REVISION of this file.

    Parameters
    ----------
    crate : a pycrates crate
        The block to use.

    Returns
    -------
    revision : int
        The revision number.

    Notes
    -----
    At the moment the revision number is extracted from the file
    name (the Nxxx part) rather than using the REVISION keyword as
    it is not clear to me the latter has much discriminatory power.
    The file name (not including path) must begin with the STACK_ID
    value and then followed with "Nddd_" (where d represents '0'
    to '9').

    At the time of writing 2018-06-19) Joe has confirmed that mrgsrc3
    files are not guaranteed to have the correct REVISION value,
    although this should be fixed once the stack has passed through
    the source properties pipeline.
    """

    infile = crate.get_filename()
    basename = os.path.basename(infile)

    stackname = crate.get_key_value('STACK_ID')
    assert stackname is not None, infile
    assert basename.startswith(stackname + 'N'), \
        'infile={} stackname={}'.format(infile, stackname)

    s = len(stackname) + 1
    assert basename[s + 3] == '_', \
        'infile={}  char={}'.format(infile, basename[s + 3])

    vstr = basename[s: s + 3]
    try:
        return int(vstr)
    except ValueError:
        raise ValueError("Unexpected version number " +
                         "{} in {}".format(vstr, infile))


def report_revision_difference(infile, revnum):
    """Provide a warning if infile does not match the expected revision.

    The warning is to stderr.

    Parameters
    ----------
    infile : str
        The file name to check. It must have a STACK_ID in its
        most-interesting block that matches the start of its name.
    revnum : int
        The revision number.

    """

    cr = pycrates.read_file(infile)
    foundrev = get_revision(cr)
    if foundrev == revnum:
        return
    elif foundrev > revnum:
        sys.stderr.write("WARNING: expected revision " +
                         "{} but found newer {}".format(revnum, foundrev))
    else:
        sys.stderr.write("WARNING: using an *older* " +
                         "revision {} than expected {}".format(foundrev, revnum))


def _read_hulls_from_mrgsrc3(mrgsrc3):
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

    Notes
    -----
    This used to be for external code, but has since been moved into
    read_master_hulls.
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
             cr.LIKELIHOOD.values,
             cr.POS.values,
             cr.get_column('EQSRC').values,
             cr.MAN_CODE.values.sum(axis=1)
             )

    out = []
    for cpt, eband, lhood, pos, eqsrc, mancode in zs:

        pos = validate_polygon(pos, report=False)
        eqsrc = validate_polygon(eqsrc, report=False)
        assert pos.shape == eqsrc.shape

        out.append({'stack': stack,
                    'component': cpt,
                    # I think there may still be some upper case band values
                    'band': eband.lower(),
                    'likelihood': lhood,
                    'mancode': mancode > 0,
                    'pos': pos.copy(),
                    'eqpos': eqsrc.copy(),
                    'transform': tr,
                    })

    return out


def read_master_hulls(chsfile, mrgsrc3dir):
    """Read in hull information.

    Parameters
    ----------
    chsfile : str
        The name of the file created by chs_create_initial_masters.py
        (the master hulls and their components for an ensemble).
    mrgsrc3dir : str
        The location of the directory containing the mrgsrc3 file
        for this ensemble. The version is expected to match that used
        to create the chsfile; a warning is displayed if it is not.

    Returns
    -------
    hullmatch, hulllist, metadata : dict, dict, dict
        The contents of the HULLMATCH and HULLLIST block, and metadata
        about the file (e.g. ensemble and CHSVER value). The
        hullmatch and hulllist dicts have keys of Master_Id.
        The hullmatch contains information on the stack-level CHS
        polygons and transform.

    Notes
    -----
    The component values *are corrected* for the COMPZERO value
    (if set, if not set it is taken to be 0), and the COMPZERO
    value is included in the output for each hull.
    """

    ds = pycrates.CrateDataset(chsfile, mode='r')
    cr = ds.get_crate('HULLMATCH')

    hullmatch = {}

    # Map from stack id to the "number" of the stack in the ensemble
    ensemblemap = {}
    for i in range(cr.get_key_value('STKIDNUM')):
        key = 'STKID{:03d}'.format(i)
        stkid = cr.get_key_value(key)
        assert stkid not in ensemblemap, stkid
        ensemblemap[stkid] = i

    # Need to support old data during testing
    compzero = cr.get_key_value('COMPZERO')
    if compzero is None:
        compzero = 0

    # As the same stack can appear multiple times, store the
    # data from the mrgsrc3 files in a dictionary. The assumption
    # is that this is not going to eat up too much memory.
    #
    mrgsrc3files = {}

    zs = zip(cr.Master_Id.values,
             cr.STACKID.values,
             cr.COMPONENT.values,
             cr.EBAND.values,
             cr.LIKELIHOOD.values,
             cr.MAN_CODE.values,
             cr.MRG3REV.values,
             cr.INCLUDE_IN_CENTROID.values,
             cr.STKSVDQA.values,
             cr.Match_Type.values)
    for vals in zs:
        mid, stackid, component, eband, lhood, man_code, \
            mrg3rev, incl_cen, svdqa, mtype = vals

        # Check there's no upper/lower-case confusion here.
        assert eband in "busmhw", eband

        # man_code is a single-element array; convert to a scalar
        # (and add a check in case this ever changes)
        #
        assert man_code.shape == (1,)
        man_code = man_code[0]

        try:
            mrgsrc3data = mrgsrc3files[stackid]
        except KeyError:
            filename = find_mrgsrc3(stackid, mrgsrc3dir)
            report_revision_difference(filename, mrg3rev)

            mrgsrc3files[stackid] = _read_hulls_from_mrgsrc3(filename)
            mrgsrc3data = mrgsrc3files[stackid]

        matches = [m3 for m3 in mrgsrc3data
                   if m3['component'] == component]
        assert len(matches) == 1
        m3 = matches[0]

        # Could check for differences between the mhull file and
        # the mrgsrc3 file, but already have a version check (which
        # does not mean the values are different), so worry about
        # adding such a check if it becomes useful.
        #
        store = {'master_id': mid,
                 'stack': stackid,
                 'component': component - compzero,
                 'compzero': compzero,
                 'eband': eband,
                 'likelihood': lhood,

                 # man_code is the actual code, mancode is a flag
                 # indicating whether 0 or not (not really needed but
                 # left in from amalgamating code)
                 'man_code': man_code,
                 'mancode': man_code > 0,

                 'mrg3rev': mrg3rev,
                 'include_in_centroid': incl_cen,
                 'stksvdqa': svdqa,
                 'match_type': mtype,

                 # polygon data
                 'pos': m3['pos'],
                 'eqpos': m3['eqpos'],
                 'transform': m3['transform']
                 }

        try:
            hullmatch[mid].append(store)
        except KeyError:
            hullmatch[mid] = [store]

    cr = ds.get_crate('HULLLIST')

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


def polygon_centroid(xs, ys):
    """Return the centroid of the polygon.

    Parameters
    ----------
    xs, ys : NumPy array
        The vertexes of the polygon, which is assumed to be closed
        (i.e. xs[0] == xs[-1], ys[0] == ys[-1]). The arrays must
        have the same length and contain at least 4 points.

    Returns
    -------
    x0, y0 : float
        The centroid of the polygon, calculated using [1]_.

    References
    ----------

    .. [1] https://en.wikipedia.org/wiki/Centroid#Centroid_of_a_polygon

    """

    # use j to mean i + 1
    #
    xi = xs[:-1]
    xj = xs[1:]
    yi = ys[:-1]
    yj = ys[1:]

    # Note that I am not simplifying the constants here as the
    # reduction in exection time should not be significant.
    #
    term = xi * yj - xj * yi
    area = term.sum() / 2.0
    denom = 6 * area

    cx = ((xi + xj) * term).sum() / denom
    cy = ((yi + yj) * term).sum() / denom

    return cx, cy
