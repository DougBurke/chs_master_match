"""
Utility routines.
"""

import glob
import json
import os
import sys
import time

import six

from collections import defaultdict

# import six

import numpy as np

import pycrates
import region

from coords.format import deg2ra, deg2dec

import chs_status


def have_directory_write_access(dirname):
    """Do we have write access to the given directory?

    Parameters
    ----------
    dirname: str
        The name of the directory to check.

    Returns
    -------
    flag : bool

    Notes:
    ------
    Based on https://stackoverflow.com/a/2113511
    """

    return os.access(dirname, os.W_OK | os.X_OK)


def log(msg, level='LOG'):
    """Log a message.

    Parameters
    ----------
    msg: str
        The message
    """

    print("{}: {}".format(level, msg))


def dbg(msg):
    """Log a debug message.

    Parameters
    ----------
    msg: str
        The message
    """

    log(msg, level='DEBUG')


def warn(msg):
    """Log a warning message.

    Parameters
    ----------
    msg: str
        The message
    """

    log(msg, level='WARNING')


def errlog(msg):
    """Log an error message.

    Parameters
    ----------
    msg: str
        The message
    """

    log(msg, level='ERROR')


_logonce = set([])


def logonce(msg):
    """Display the message to STDOUT if it hasnt' already been reported.

    The message is preceeded by 'LOG: '.

    Parameters
    ----------
    msg : str
        The message to display (actually can be anything that can
        be converted to a string and has an equality test).

    """

    if msg in _logonce:
        return

    log(msg)
    _logonce.add(msg)


def touch_file(filename):
    """Morally similar to the touch command.

    Based on https://stackoverflow.com/a/1160227 since I am
    still supporting Python 2.7 and do not want to add another
    OTS package ('pip install pathlib').

    Parameters
    ----------
    filename : str
        The path to the file to "touch".
    """

    with open(filename, 'a'):
        os.utime(filename, None)


def find_single_match(pat):
    """Return the file name matching the pattern.

    If there's a single match, return it, otherwise error (so
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


def make_terminal_name(ensemble):
    """The file name used to indicate that the ensemble is complete.

    Parameters
    ----------
    ensemble : str
        The ensemble name.

    Returns
    -------
    filename : str
        The file name (no path is included).

    """

    return "COMPLETED.{}".format(ensemble)


def make_mhull_name(ensemble, revision=None):
    """What is the name of the master hull file?

    If the revision is not given then a wild card is added ("*").

    Parameters
    ----------
    ensemble : str
        The ensemble name.
    revision : int or None, optional
        The revision number (if known).

    Returns
    -------
    filename : str
        The name of the mhull file (with no path).
    """

    if revision is None:
        revstr = "*"
    else:
        revstr = "{:03d}".format(int(revision))

    return 'master_hulls.{}.v{}.fits'.format(ensemble, revstr)


def find_mhulls(indir, ensemble):
    """Return mhull files for this ensemble, indexed by CHSVER.

    Note that the input directory is not the usual "datadir" path
    as it includes the ensemble name (for the standard data
    arrangement).

    Parameters
    ----------
    indir : str
        The path to an ensemble directory, containing mhull files
        for the given ensemble.

    Returns
    -------
    ans : [(int, str)]
        The version (CHSVER) and file name for all the mhull files
        for this ensemble. There must be at least 1. The return
        list is sorted so that the versions are in numerical order
        (increasing); i.e. ans[0][0] = 1, ans[-1] is the latest
        version.

    See Also
    --------
    find_field_jsons, make_mhull_name

    """

    pat = os.path.join(indir,
                       make_mhull_name(ensemble))
    matches = glob.glob(pat)
    if len(matches) == 0:
        raise IOError("No mhull file found matching {}".format(pat))

    vers = []
    for match in matches:
        ver = pycrates.read_file(match).get_key_value('CHSVER')
        if ver is None:
            raise IOError("Missing CHSVER in {}".format(match))

        vers.append(int(ver))

    nver = len(vers)
    ncheck = len(set(vers))
    if nver != ncheck:
        raise IOError("Versions not unique: {}".format(vers))

    out = sorted(list(zip(vers, matches)), key=lambda x: x[0])
    return out


def find_field_jsons(indir, ensemble):
    """Return JSON field files for this ensemble, indexed by CHSVER.

    Note that the input directory is not the usual "userdir" path
    as it includes the ensemble name.

    Parameters
    ----------
    indir : str
        The path to an ensemble directory, containing JSON files
        for the given ensemble created by the server.

    Returns
    -------
    ans : [(int, str)]
        The version (CHSVER) and file name for all the field JSON
        files for this ensemble. There must be at least 1. The return
        list is sorted so that the versions are in numerical order
        (increasing); i.e. ans[0][0] = 1, ans[-1] is the latest
        version.

    See Also
    --------
    find_mhulls, make_field_name_json

    Notes
    -----
    At present this is unused.

    """

    pat = os.path.join(indir,
                       make_field_name_json(ensemble))
    matches = glob.glob(pat)
    if len(matches) == 0:
        raise IOError("No mhull file found matching {}".format(pat))

    vers = []
    for match in matches:

        jcts = read_json(match)
        ver = int(jcts['revision'])
        vers.append(int(ver))

    nver = len(vers)
    ncheck = len(set(vers))
    if nver != ncheck:
        raise IOError("Versions not unique: {}".format(vers))

    out = sorted(list(zip(vers, matches)), key=lambda x: x[0])
    return out


def make_field_name_json(ensemble, revision=None):
    """What is the JSON file storing the field data?

    If the revision is not given then a wild card is added ("*").

    Parameters
    ----------
    ensemble : str
        The ensemble name.
    revision : int or None, optional
        The revision number (if known)

    Returns
    -------
    filename : str
        The name of the JSON file (with no path).
    """

    if revision is None:
        revstr = "*"
    else:
        revstr = "{:03d}".format(int(revision))

    return 'field.{}.v{}.json'.format(ensemble, revstr)


def make_hull_name_json(ensemble, masterid, revision):
    """What is the JSON file storing the master hull data?

    Parameters
    ----------
    ensemble : str
        The ensemble name.
    masterid : int or None.
        The master id value.
    revision : int
        The revision number.

    Returns
    -------
    filename : str
        The name of the JSON file (with no path). If masterid
        is None then a wild card is used.
    """

    # assert isinstance(revision, six.string_types)
    revision = int(revision)

    if masterid is None:
        mstr = "*"
    else:
        mstr = "{:03d}".format(masterid)

    return 'hull.{}.{}.v{:03d}.json'.format(ensemble,
                                            mstr,
                                            revision)


def make_poly_name_json(ensemble, masterid, revision):
    """What is the JSON file storing the master hull polygon?

    Parameters
    ----------
    ensemble : str
        The ensemble name.
    masterid : int
        The master id value.
    revision : int
        The revision number.

    Returns
    -------
    filename : str
        The name of the JSON file (with no path).
    """

    # assert isinstance(revision, six.string_types)
    revision = int(revision)

    return 'poly.{}.{:03d}.v{:03d}.json'.format(ensemble,
                                                masterid,
                                                revision)


def make_component_name_json(ensemble, stack, component, revision):
    """What is the JSON file storing the stack-level data?

    Parameters
    ----------
    ensemble : str
        The ensemble name.
    stack : str
        The stack name.
    component : int
        The component number.
    revision : int
        The revision number.

    Returns
    -------
    filename : str
        The name of the JSON file (with no path).
    """

    # assert isinstance(revision, six.string_types)
    revision = int(revision)

    return 'cpt.{}.{}.{:02d}.v{:03d}.json'.format(ensemble,
                                                  stack,
                                                  component,
                                                  revision)


def make_component_region_name(stack, component, revision):
    """What is the reg file for this component?

    Parameters
    ----------
    stack : str
        The stack name.
    component : int
        The component number.
    revision : int
        The revision number.

    Returns
    -------
    filename : str
        The name of the region file (with no path).
    """

    # assert isinstance(revision, six.string_types)
    revision = int(revision)

    # Note: no zero-padding for the component number
    return 'stack.{}.{}.v{:03d}.reg'.format(stack,
                                            component,
                                            revision)


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


def convert_to_key(ensemblemap, stack, component):
    """Convert stack name and component to a "key".

    This key can be useful when we only want to index by a single
    value, rather than a tuple (e.g. passing information via
    JSON/Javascript).

    Parameters
    ----------
    ensemblemap : dict
        The keys are stack names, the values are the integer
        value for the stack in the ensemble (i.e. the XXX value
        from the STKIDXXX keyword).

    Returns
    -------
    key : str
        The key - at present a 3-character integer, a '.', and then
        a 2-character integer, where the numbers are 0-padded to the
        left.
    """

    try:
        return '{:03d}.{:02d}'.format(ensemblemap[stack],
                                      component)
    except KeyError:
        raise ValueError("Unrecognised stack {}".format(stack) +
                         " in:\n{}".format(ensemblemap))


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

    hullmatch = defaultdict(list)

    # Map from stack id to the "number" of the stack in the ensemble
    ensemblemap = {}
    stacklist = []
    for i in range(cr.get_key_value('STKIDNUM')):
        key = 'STKID{:03d}'.format(i)
        stkid = cr.get_key_value(key)
        assert stkid not in ensemblemap, stkid
        stacklist.append(stkid)
        ensemblemap[stkid] = i

    compzero = cr.get_key_value('COMPZERO')
    assert compzero is not None

    svdqafile = cr.get_key_value('SVDQAFIL')
    centroidfile = cr.get_key_value('CENFILE')
    assert svdqafile is not None
    assert centroidfile is not None

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
        # In subsequent versions it can be back to a scalar, so
        # support this. This is because it gets written out as
        # Int4 rather than Byte.
        #
        if man_code.shape == ():
            pass
        elif man_code.shape == (1,):
            man_code = man_code[0]
        else:
            raise ValueError("man_code.shape = {}".format(man_code.shape))

        try:
            mrgsrc3data = mrgsrc3files[stackid]
        except KeyError:
            filename = find_mrgsrc3(stackid, mrgsrc3dir)
            report_revision_difference(filename, mrg3rev)

            mrgsrc3files[stackid] = _read_hulls_from_mrgsrc3(filename)
            mrgsrc3data = mrgsrc3files[stackid]

        key = convert_to_key(ensemblemap, stackid, component)

        matches = [m3 for m3 in mrgsrc3data
                   if m3['component'] == component]
        assert len(matches) == 1
        m3 = matches[0]

        # Could check for differences between the mhull file and
        # the mrgsrc3 file, but already have a version check (which
        # does not mean the values are different), so worry about
        # adding such a check if it becomes useful.
        #
        # Explicitly convert from NumPy to Python types for some
        # integers, since np.intXX types tend not to be JSON
        # serializable.
        #
        store = {'master_id': int(mid),
                 'stack': stackid,
                 'component': int(component - compzero),
                 'key': key,
                 'compzero': int(compzero),
                 'eband': eband,
                 'likelihood': lhood,

                 # man_code is the actual code, mancode is a flag
                 # indicating whether 0 or not (not really needed but
                 # left in from amalgamating code)
                 #
                 'man_code': int(man_code),
                 'mancode': man_code > 0,

                 'mrg3rev': int(mrg3rev),
                 'include_in_centroid': incl_cen,
                 'stksvdqa': svdqa,
                 'match_type': mtype,

                 # polygon data
                 'pos': m3['pos'],
                 'eqpos': m3['eqpos'],
                 'transform': m3['transform']
                 }

        hullmatch[mid].append(store)

    cr = ds.get_crate('HULLLIST')

    hulllist = {}
    zs = zip(cr.Master_Id.values,
             cr.STATUS.values,
             cr.BASE_STK.values,
             cr.NVERTEX.values,
             cr.EQPOS.values)
    for mid, status, base_stk, nvertex, eqpos in zs:
        assert mid not in hulllist

        hulllist[mid] = {'master_id': int(mid),
                         'status': status,
                         'base_stk': base_stk,
                         'eqpos': eqpos[:, :nvertex]}

    metadata = {'ensemble': cr.get_key_value('ENSEMBLE'),
                'ensemblemap': ensemblemap,
                'stacks': stacklist,
                'svdqafile': svdqafile,
                'centroidfile': centroidfile,
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


# JSON code

def read_json(infile):
    """Return data stored as JSON in the input file.

    Parameters
    ----------
    infile : str
        The full path to the JSON file.

    Returns
    -------
    ans : dict or array or None
        The contents, or None if there was an error.
    """

    try:
        with open(infile, 'r') as fh:
            cts = fh.read()
            jcts = json.loads(cts)

    except Exception as exc:
        errlog("error reading JSON from {}\n{}".format(infile, exc))
        return None

    return jcts


def setup_user_setting(store, key, stringval=True):
    """Add in user-level version.

    This is highly-specialized. It changes the key value to
    be a dictionary with 'proposed' and 'user' keywords,
    storing the current value under 'proposed' and setting
    'user' to None.

    Parameters
    ----------
    store : dict
        A dictionary which is assumed to contain the supplied key.
    key : dict_key
        The key value.
    stringval : bool, optional
        If True then the value must be a string.

    Returns
    -------
    flag : bool
        If False then the stored value was not a string so nothing
        has been changed. It is expected that this will cause
        down stream to trigger an error handler. This is only
        set to False if stringval is True.

    Notes
    -----
    If the key does not exist in the input store then it its "proposed"
    value is set to '' (stringval is True) or None (otherwise).

    TODO: validate that setting value to '' is actually what we want,
    since there is some JSON code that expects it to be None/null.
    """

    try:
        v = store[key]
    except KeyError:
        warn("NO {} field".format(key))
        if stringval:
            v = ''
        else:
            v = None

    if stringval and not isinstance(v, six.string_types):
        errlog("{} is not a string but {}".format(key, v))
        return False

    store[key] = {'proposed': v, 'user': None}
    return True


def get_user_setting(store, key):
    """Return the setting for this key.

    This is a key that has been through the setup_user_setting
    routine, so has a 'proposed' and 'user' variant. The
    'user' setting is returned unless it is None or '', in which
    case the 'proposed' setting is used (and a check is made to
    ensure it is not None, but '' *is* allowed).

    Parameters
    ----------
    store : dict
        A dictionary which is assumed to contain the supplied key.
    key : dict_key
        The key value.

    Returns
    -------
    value
        The value for this key.

    Notes
    -----
    The code base has not been clear about what a value='' means;
    so it is not clear what should be done if the proposed value
    is ''. For now this is being allowed.
    """

    try:
        v = store[key]
    except KeyError:
        raise KeyError("Store does not contain key={}".format(key))

    # Always extract both as a check it is a user-settable value.
    # Do not bother trying to make a "user friendly" error message.
    #
    proposed = v['proposed']
    user = v['user']

    if user is not None and user != '':
        return user

    # if proposed is None or proposed == '':
    if proposed is None:
        raise ValueError("Unexpected user/proposed values for " +
                         "key={} in store=\n{}".format(key, store))

    return proposed


def read_ensemble_hull_json(datadir, userdir,
                            ensemble, mid, revision):
    """Return JSON-stored data for the masters.

    Parameters
    ----------
    datadir : str
        The path to the directory containing the ensemble-level
        products.
    userdir : str
        The path to the directory containing the user's decisions.
    ensemble : str
        The ensemble.
    mid : int
        The master id.
    revision : str
        Expected to be 001, ...

    Returns
    -------
    ans : dict or None
        The contents, or None if there was an error.

    """

    filename = make_hull_name_json(ensemble, mid, revision)

    infile = os.path.join(datadir, ensemble, filename)
    jcts = read_json(infile)
    if jcts is None:
        return None

    # Setup for user information.
    #
    userkeys = ['usernotes', 'useraction', 'lastmodified']
    for key in userkeys:
        if not setup_user_setting(jcts, key):
            return None

    # Now add in any user information
    #
    infile = os.path.join(userdir, ensemble, filename)

    # This is an optional file, so avoid warning messages in the log
    # if we can help it.
    if not os.path.exists(infile):
        return jcts

    ucts = read_json(infile)
    if ucts is None:
        return jcts

    # Only override with the user values
    #
    for key in userkeys:
        if key in ucts:
            jcts[key]['user'] = ucts[key]

    return jcts


def read_ensemble_json(datadir, userdir, ensemble):
    """Return JSON data from the ensemble.

    The returned structure contains all the versions for this
    ensemble. The "current" version is taken to be the highest
    version number.

    Parameters
    ----------
    datadir : str
        The path to the directory containing the ensemble-level
        products. These are the "proposed" products.
    userdir : str
        The path to the directory containing the user's decisions.
    ensemble : str
        The ensemble.

    Returns
    -------
    ans : dict or None
        The contents, or None if there was an error.

    Notes
    -----

    state['latest_version'] => str (001, 002, ..)
    state['versions'] => dict, keys are version
          ['nstacks']
          ['nmasters']  == len(masters)
          ['masters'] = list of dict
                  ['masterid'] => str (001, ...)
                  ['ncpts'] == len(cpts)
                  ---- ['cpts'] => list of ?
                  ['usernotes'] =>
                  ['useraction'] =>
                  ['lastmodified'] =>

           ['usernotes']
           ['lastmodified']

    The usernotes, useraction, lastmodified, and status fields are
    dicts with two keys:
      proposed
      user
    which allows the user to "over ride" the proposed value.

hmmm, the JSON data has {"status": "todo", "stackmap": {"acisfJ1705367m403832_001": 3, "acisfJ1704041m414416_001": 1, "acisfJ1702545m412821_001": 0, "acisfJ1705559m410515_001": 4, "acisfJ1704448m410953_001": 2}, "nmasters": 1, "name": "ens0000900_001", "lastmodified": "", "ncpts": 2, "usernotes": "", "nstacks": 2, "revision": "001"}


    """

    # Process the proposed settings first
    #
    pat = "field.{}.*.json".format(ensemble)
    inpat = os.path.join(datadir, ensemble, pat)
    matches = glob.glob(inpat)
    if len(matches) == 0:
        errlog("no field.json files found for ensemble " +
               "{} - {}".format(ensemble, inpat))
        return None

    store = {'versions': {}}
    for match in matches:
        jcts = read_json(match)
        if jcts is None:
            continue

        try:
            v = jcts['revision']
        except KeyError:
            log("missing revision keyword in {}".format(match))
            continue

        # this should not happen, so do not worry too much about the
        # error handler
        if v in store['versions']:
            log("multiple revision={} in {}".format(v, match))
            continue

        # TODO: the following loop requires that the masters are always
        #       numbered 1, ..., n (with no gaps).
        #       Is this always guaranteed? Probably not.
        #
        hulls = []
        for mid in range(1, jcts['nmasters'] + 1):
            hull = read_ensemble_hull_json(datadir, userdir,
                                           ensemble, mid, v)
            if hull is None:
                # if there is a problem reading in a single hull, then
                # bail out for the whole thing
                return None

            hulls.append(hull)

        if len(hulls) == 0:
            log("no masters for revision " +
                "{} in {}".format(v, match))
            continue

        if 'masters' in jcts:
            log("overwriting masters setting in " +
                "{}".format(match))

        jcts['masters'] = hulls
        store['versions'][v] = jcts

        # Override those keys that need proposed/user versions.
        #
        for key in ['usernotes', 'lastmodified', 'status']:
            if not setup_user_setting(jcts, key):
                return None

    revs = list(store['versions'].keys())
    if len(revs) == 0:
        errlog("no JSON data read from ensemble: " +
               "{} {}".format(datadir, ensemble))
        return None

    revs = sorted(revs, key=int, reverse=True)
    store['latest_version'] = revs[0]

    # Now check for user overrides: at present only at the
    # field level.
    #
    inpat = os.path.join(userdir, ensemble, pat)
    matches = glob.glob(inpat)
    if len(matches) == 0:
        return store

    # Assume very-limited metadata here.
    #
    for match in matches:
        jcts = read_json(match)

        revision = jcts['revision']

        try:
            base = store['versions'][revision]
        except KeyError:
            warn("{} has invalid version {}".format(match,
                                                    revision))
            continue

        for k in ['lastmodified', 'usernotes', 'status']:
            try:
                base[k]['user'] = jcts[k]
            except KeyError:
                # for now do not require the user to have all fields
                # set.
                pass

    return store


def read_ensemble_status(datadir, userdir, ensemble, revision):
    """What is the ensemble status?

    What is the current status of the ensemble? This is based on
    read_ensemble_json.

    Parameters
    ----------
    datadir : str
        The path to the directory containing the ensemble-level
        products. These are the "proposed" products.
    userdir : str
        The path to the directory containing the user's decisions.
    ensemble : str
        The ensemble.
    revision : str
        The revision value (in 3-digit form).

    Returns
    -------
    status : {'todo', 'review', 'completed', 'delete', 'done', 'unknown'}

    """

    # Try the user and then proposed settings.
    #
    pat = make_field_name_json(ensemble, revision)

    # pick an out-of-bounds value; normally I'd use None but
    # this is a possible value, so use a numeric value.
    #
    not_present = -1

    def lookin(indir):

        infile = os.path.join(indir, ensemble, pat)
        if not os.path.exists(infile):
            return not_present

        try:
            jcts = read_json(infile)
        except IOError:
            return not_present

        try:
            return jcts['status']
        except KeyError:
            return not_present

    ans = lookin(userdir)
    if ans is not None and ans != not_present:
        return ans

    ans = lookin(datadir)
    if ans is not None and ans != not_present:
        return ans

    # The following logic is not clear
    if ans is None:
        warn("status=None for ensemble " +
             "{} version {}".format(ensemble, revision))
        return chs_status.UNKNOWN

    errlog("no status for ensemble " +
           "{} version {}".format(ensemble, revision))
    return chs_status.UNKNOWN


def read_component_hull_json(datadir, userdir, ensemble,
                             stack, component, revision):
    """What do we have stored for this component hull?

    Parameters
    ----------
    datadir : str
        The path to the directory containing the ensemble-level
        products. These are the "proposed" products.
    userdir : str
        The path to the directory containing the user's decisions.
    ensemble : str
        The ensemble.
    stack : str
        The stack identifier (e.g. 'acisf...' or 'hrcf...')
    component : int
        The component number of the stack-level hull.
    revision : str
        The revision value (in 3-digit form).

    Returns
    -------
    ans : dict or None
        The contents, or None if there was an error.

    """

    filename = make_component_name_json(ensemble,
                                        stack,
                                        component,
                                        revision)

    infile = os.path.join(datadir, ensemble, filename)
    jcts = read_json(infile)
    if jcts is None:
        return None

    # Setup for user information.
    #
    usermodkeys = ['master_id', 'include_in_centroid', 'lastmodified']
    usermodstrings = [False, False, True]

    for key, flag in zip(usermodkeys, usermodstrings):
        if not setup_user_setting(jcts, key, stringval=flag):
            return None

    # Now add in any user information
    #
    infile = os.path.join(userdir, ensemble, filename)

    # This is an optional file, so avoid warning messages in the log
    # if we can help it.
    if not os.path.exists(infile):
        return jcts

    ucts = read_json(infile)
    if ucts is None:
        return jcts

    for key in usermodkeys:
        if key in ucts:
            jcts[key]['user'] = ucts[key]

    return jcts


def save_datatable(userdir, data):
    """Save the data table details.

    This could be included in the summary page, but for now
    separate it out.

    Parameters
    ----------
    userdir : str
        The location for the user-stored data.
    data : dict
        The JSON dictionary containing the elements to write out.
    """

    outname = 'datatable.json'
    outfile = os.path.join(userdir, outname)

    with open(outfile, 'w') as fh:
        fh.write(json.dumps(data))


def save_summary(userdir, data):
    """Save the summary details.

    Parameters
    ----------
    userdir : str
        The location for the user-stored data.
    data : dict
        The JSON dictionary containing the elements to write out.
    """

    outname = 'summary.json'
    outfile = os.path.join(userdir, outname)

    store = {"lastmodified": time.asctime(),
             "usernotes": data['usernotes']}
    with open(outfile, 'w') as fh:
        fh.write(json.dumps(store))


def save_ensemble(userdir, data):
    """Save the ensemble-level details.

    It will not over-write an existing file.

    Parameters
    ----------
    userdir : str
        The location for the user-stored data.
    data : dict
        The JSON dictionary containing the elements to write out.
        Only a restricted subset is written out.

    Returns
    -------
    outfile : str
        The name of the file that was created.
    """

    ensemble = data['name']
    version = data['revision']  # this is in string form, 0 padded

    outdir = os.path.join(userdir, ensemble)
    outname = make_field_name_json(ensemble, version)
    outfile = os.path.join(outdir, outname)
    if os.path.exists(outfile):
        raise IOError("Output file already exists: {}".format(outfile))

    if os.path.exists(outdir):
        if not os.path.isdir(outdir):
            raise IOError("Exists but not a directory! {}".format(outdir))
    else:
        os.mkdir(outdir)

    store = {"name": ensemble,
             "lastmodified": time.asctime(),
             "usernotes": data['usernotes'],
             "status": data['status'],
             "revision": version}

    for key in ['stackmap', 'nmasters', 'nstacks']:
        if key in data:
            store[key] = data[key]

    with open(outfile, 'w') as fh:
        fh.write(json.dumps(store))

    return outfile


def save_master(userdir, data):
    """Save the master-level details.

    It will not over-write an existing file. It is used for both the
    proposed and user versions of the file, and at present not all
    information is required for the user version (since it is meant
    to overwrite some, but not all, of the proposed values).

    Parameters
    ----------
    userdir : str
        The location for the user-stored data.
    data : dict
        The JSON dictionary containing the elements to write out.
        Only a restricted subset is written out.

    Returns
    -------
    outfile : str
        The name of the file that was created.
    """

    ensemble = data['ensemble']
    version = data['revision']  # this is in string form, 0 padded
    masterid = int(data['masterid'])

    outdir = os.path.join(userdir, ensemble)
    outname = make_hull_name_json(ensemble, masterid, version)
    outfile = os.path.join(outdir, outname)
    if os.path.exists(outfile):
        raise IOError("Output file already exists: {}".format(outfile))

    if os.path.exists(outdir):
        if not os.path.isdir(outdir):
            raise IOError("Exists but not a directory! {}".format(outdir))
    else:
        os.mkdir(outdir)

    store = {"ensemble": ensemble,
             "masterid": "{:03d}".format(masterid), # is this sensible?
             "usernotes": data['usernotes'],
             "useraction": data['useraction'],
             "revision": version}

    if 'lastmodified' in data:
        store['lastmodified'] = data['lastmodified']
    else:
        store['lastmodified'] = time.asctime()

    # do we ever want nstacks?
    for key in ['ncpts', 'nstacks', 'npages', 'usernotes']:
        if key in data:
            store[key] = data[key]

    # hard code npages if not set
    #
    if 'ncpts' in store and 'npages' not in store:
        pagesize = 9
        ncpts = store['ncpts']
        npages = ncpts // pagesize
        if ncpts % pagesize > 0:
            npages += 1

        store['npages'] = npages

    with open(outfile, 'w') as fh:
        fh.write(json.dumps(store))

    return outfile


def save_master_poly(userdir, data):
    """Save the masterhull polygon(s).

    Parameters
    ----------
    userdir : str
        The location for the user-stored data.
    data : dict
        The JSON dictionary containing the elements to write out.
    """

    ensemble = data['ensemble']
    version = data['revision']  # this is in string form, 0 padded
    masterid = data['masterid']

    outdir = os.path.join(userdir, ensemble)
    if os.path.exists(outdir):
        if not os.path.isdir(outdir):
            raise IOError("Exists but not a directory! {}".format(outdir))
    else:
        os.mkdir(outdir)

    outname = make_poly_name_json(ensemble, masterid, version)
    outfile = os.path.join(outdir, outname)

    store = {"ensemble": ensemble,
             "masterid": masterid,
             "lastmodified": time.asctime(),
             "polygons": data['polygons'],
             "revision": version}
    with open(outfile, 'w') as fh:
        fh.write(json.dumps(store))


def save_component(userdir, data):
    """Save the stack-level hull details.

    It will not over-write an existing file.

    Parameters
    ----------
    userdir : str
        The location for the user-stored data.
    data : dict
        The JSON dictionary containing the elements to write out.

    Returns
    -------
    outfile : str
        The name of the file that was created.
    """

    ensemble = data['ensemble']
    version = data['revision']  # this is in string form, 0 padded
    stack = data['stack']
    component = data['component']

    outdir = os.path.join(userdir, ensemble)
    outname = make_component_name_json(ensemble,
                                       stack,
                                       component,
                                       version)
    outfile = os.path.join(outdir, outname)
    if os.path.exists(outfile):
        raise IOError("Output file already exists: {}".format(outfile))

    if os.path.exists(outdir):
        if not os.path.isdir(outdir):
            raise IOError("Exists but not a directory! {}".format(outdir))
    else:
        os.mkdir(outdir)

    # Could just leave this to copying over all the elements in data,
    # but at least this way there's some documentation of important
    # and required fields (could be checked in other ways though).
    #
    store = {"ensemble": ensemble,
             "stack": stack,
             "component": component,
             "lastmodified": time.asctime(),
             "revision": version,
             # DO WE NEED TO GET THE USER SETTING?
             "master_id": data['master_id'],
             "include_in_centroid": data['include_in_centroid']
             }

    # Copy over the other elements
    #
    for k, v in data.items():
        if k in store:
            continue

        if isinstance(v, dict):
            # this repeats some of the work already done
            v = get_user_setting(data, k)

        store[k] = v

    with open(outfile, 'w') as fh:
        fh.write(json.dumps(store))

    return outfile


# Crates helpers

def add_header(cr, hdrvals):
    """Add the header keywords to the crate.

    Parameters
    ----------
    cr : pycrates.Crate
    hdrvals : sequence of (name, value, desc) triples
    """

    for name, value, desc in hdrvals:
        key = pycrates.CrateKey()
        key.name = name
        key.value = value
        key.desc = desc
        cr.add_key(key)


# Should send in the time value
_header = {}


def add_standard_header(cr, creator=None, revision=1):
    """Add standard header keywords.

    The DATE value is set the first time that this routine is called,
    whatever file is being created. This may or may not be
    rather silly.
    """

    if 'timestr' not in _header:
        _header['timestr'] = time.strftime("%Y-%m-%dT%H:%M:%S")

    timestr = _header['timestr']

    hdrvals = []
    if creator is not None:
        hdrvals.append(('CREATOR', creator,
                        'tool that created this output'))

    hdrvals.extend([('DATE', timestr,
                     'Date and time of file creation'),
                    ('CHSVER', revision,
                     'The revision number for the file')
                    ])

    add_header(cr, hdrvals)


def add_col(cr, name, values,
            desc=None,
            unit=None):
    """Add a column to the crate.

    Parameters
    ----------
    cr : pycrates.TABLECrate
    name : str
    values : array_like
    desc : str or None, optional
        The description.
    unit : str or None, optional
        The units for the column.

    """

    col = pycrates.CrateData()
    col.name = name
    col.values = values
    if desc is not None:
        col.desc = desc
    if unit is not None:
        col.unit = unit

    cr.add_column(col)


# Create the master hull file (mhull).
#
def read_svdqafile(infile):
    """Returns the set of stacks that went to SVD QA.

    Parameters
    ----------
    infile : str
        The name of the file. The first column is used, and it
        is assumed to be case-sensitive.

    Returns
    -------
    stackids : set of str
        The stacks that have been to SVD QA.

    """

    cr = pycrates.read_file(infile)
    if cr.get_nrows() < 1:
        raise IOError("No data read from: {}".format(infile))

    # Explicit conversion to a Python string
    return set([str(v) for v in cr.get_column(0).values])


def read_centroidfile(infile):
    """The set of stack-level hulls to exclude from centroid calculation.

    The file is assumed to have columns "stack", "cpt", and
    "use_cen". Only those with "use_cen" set to False are returned
    here.

    Parameters
    ----------
    infile : str
        The name of the file.

    Returns
    -------
    hulls : set of (stackid, component)
        The stack-level hulls to exclude from centroid calculation.

    Notes
    -----
    Although technically there could be a case where no stack-level
    hulls are to be excluded, we know this is not the case here,
    so an error is raised if no such rows are found (as a sanity
    check).

    """

    cr = pycrates.read_file(infile)
    if cr.get_nrows() < 1:
        raise IOError("No data read from: {}".format(infile))

    # do filtering here rather than with a DM filter as not 100%
    # convinced this works correctly in CIAO 4.9 (there have been
    # problems with string filtering).
    #
    out = set([])
    for stack, cpt, flag in zip(cr.stack.values,
                                cr.cpt.values,
                                cr.use_cen.values):
        if flag != "False":
            continue

        key = stack, cpt
        out.add(key)

    if len(out) == 0:
        raise IOError("No excluded data found in {}".format(infile))

    return out


def create_mhull_file(ensemble, revision, outfile,
                      hullmatch, hulllist, header,
                      creator=None):
    """Create the "CHS mst3" file.

    Parameters
    ----------
    ensemble : string
        The ensemble value, written to the header as the ENSEMBLE
        keyword.
    revision : int
        The value to write out to the header as the CHSVER
        keyword.
    outfile : string
        This file is overwritten if it exists.
    hullmatch : dict
        The contents of the HULLMATCH block. The keys are lower-case
        versions of the column names. Note that COMPONENT does not
        include the COMPZERO correction (i.e. this routine adds the
        offset).
    hulllist : dict
        The contents of the HULLLIST block. The keys are lower-case
        versions of the column names.
    header : dict
        Information used for the header (and perhaps for columns too).
        Keys are: svdqa, centroid, stacks, mstzero, and compzero.
        svdqa and centroid contain the full path to the files used
        for the STKSVDQA and INCLUDE_IN_CENTROID columns.
        stacks is a list of the stacks that form the ensemble, and to
        be written out as STKIDxxx values (in input order).
        mstzero is the value added to the Master_Id value for
        values > 0 for both HULLMATCH and HULLLIST blocks
        compzero is expected to be 0 (it is added to the COMPONENT
        values for hullmatch)
    creator : None or str, optional
        The name to use for the CREATOR field in the header.

    Notes
    -----

    It is required that the HULLMATCH block has at least one row,
    but the HULLLIST block can have 0 rows.

    """

    # As may not be set
    #
    try:
        mstzero = header['mstzero']
    except KeyError:
        mstzero = 0

    # Extract header info
    #
    extra_hdr = [('ENSEMBLE', ensemble, 'The ensemble'),
                 ('MSTZERO', mstzero,
                  'The offset applied to positive Master_Id values'),
                 ('COMPZERO', header['compzero'],
                  'The COMPONENT value for cpt=0'),
                 ('SVDQAFIL', header['svdqafile'],
                  'The stacks that went to SVD QA'),
                 ('CENFILE', header['centroidfile'],
                  'centroid input')]

    for i, stack in enumerate(header['stacks']):
        extra_hdr.append(('STKID{:03d}'.format(i),
                          stack,
                          'Member of the ensemble'))

    extra_hdr.append(('STKIDNUM', len(header['stacks']),
                      'Number of stacks in ensemble'))

    ds = pycrates.CrateDataset()

    cr = pycrates.TABLECrate()
    cr.name = 'HULLMATCH'

    add_standard_header(cr, creator=creator, revision=revision)
    add_header(cr, extra_hdr)

    # NOTE: the NHULLS column might be better in the next block,
    #       but it can be useful to know how many stack-level
    #       hulls there are in a master when looking at this data.
    #
    mids = []
    for mid in hullmatch['master_id']:
        if mid > 0:
            mid += mstzero

        mids.append(mid)

    add_col(cr, 'Master_Id', mids,
            desc='This is an internal number, do not expose')
    add_col(cr, 'NHULLS', hullmatch['nhulls'],
            desc='The number of stack-level hulls in the master')
    add_col(cr, 'STACKID', hullmatch['stackid'])

    cpts = []
    for cpt in hullmatch['component']:
        cpts.append(cpt + header['compzero'])

    add_col(cr, 'COMPONENT', cpts,
            desc='Offset by COMPZERO from MEXTSRC component value')

    add_col(cr, 'Match_Type', hullmatch['match_type'])
    add_col(cr, 'AREA', hullmatch['area'],
            unit='arcsec**2',
            desc='Area of hull excluding pixel-mask filter')
    add_col(cr, 'EBAND', hullmatch['eband'],
            desc='Energy band of hull')
    add_col(cr, 'LIKELIHOOD', hullmatch['likelihood'],
            desc='Likelihood of hull')
    add_col(cr, 'MAN_CODE', hullmatch['man_code'],
            desc='Copied from MEXTSRC block (converted to int)')
    add_col(cr, 'MRG3REV', hullmatch['mrg3rev'],
            desc='Revision of mrgsrc3 file used')
    add_col(cr, 'INCLUDE_IN_CENTROID',
            hullmatch['include_in_centroid'],
            desc='Use hull in centroid calculation?')
    add_col(cr, 'STKSVDQA', hullmatch['stksvdqa'],
            desc='Did this stack go to SVD QA?')

    ds.add_crate(cr)

    cr = pycrates.TABLECrate()
    cr.name = 'HULLLIST'

    add_standard_header(cr, creator=creator, revision=revision)
    add_header(cr, extra_hdr)

    mids = []
    for mid in hulllist['master_id']:
        if mid > 0:
            mid += mstzero

        mids.append(mid)

    add_col(cr, 'Master_Id', mids,
            desc='This is an internal number, do not expose')
    add_col(cr, 'STATUS', hulllist['status'],
            desc='Did the master-match work?')
    add_col(cr, 'BASE_STK', hulllist['base_stk'],
            desc='The stack used for SKY coord system, or NONE')

    add_col(cr, 'MANMATCH', hulllist['manmatch'],
            desc='Has the selection of stack-level hulls been changed')
    add_col(cr, 'MANREG', hulllist['manreg'],
            desc='Has the region been changed manually')

    add_col(cr, 'NVERTEX', hulllist['nvertex'],
            desc='The number of vertexes in the closed hull')

    add_col(cr, 'NSTKHULL', hulllist['nstkhull'],
            desc='The number of stack hulls that were combined')

    # NOTE: there is no POS coordinate column in this block
    #       (so can not attach a transform to it)
    #
    # Want a 3D NumPy array but assume the input is a list of
    # 2D NumPy arrays of different sizes, given by the nvertex
    # column.
    #
    nhulls = len(hulllist['nvertex'])
    if nhulls > 0:

        nmax = np.asarray(hulllist['nvertex']).max()
        eqpos = np.full((nhulls, 2, nmax), np.nan, dtype=np.float64)

        any_non_qa = False
        for i, eq in enumerate(hulllist['eqpos']):

            # Should we not include the hull polygon for status=delete?
            #
            if eq is None:
                status = hulllist['status'][i]
                assert chs_status.is_qa(status), status
                continue

            assert eq.ndim == 2, eq.shape
            assert eq.shape[0] == 2, eq.shape
            assert eq.shape[1] <= nmax, (eq.shape, nmax)

            eqpos[i, :, :hulllist['nvertex'][i]] = eq
            any_non_qa |= True

        # The nmax check is not made if all hulls are QA cases
        if any_non_qa:
            assert nmax > 2, "Expected > 2 vertex but found {}".format(nmax)

    else:
        # Can we get array with just [] here or does it have
        # to be a three-dimensional empty array? It looks like it
        # has to be nD otherwise you can get a core dump (not sure
        # what the requirements are).
        #
        eqpos = np.asarray([], dtype=np.float64).reshape(0, 2, 0)

    col = pycrates.create_vector_column('EQPOS', ['RA', 'DEC'])
    col.desc = 'The master hull vertices'
    col.unit = 'degree'
    col.values = eqpos
    cr.add_column(col)

    ds.add_crate(cr)

    # ds.write(outfile, clobber=True)
    ds.write(outfile, clobber=False)
    print("Created: {}".format(outfile))


def get_masterid_json(infile):
    """Return the value of the masterid keyword in the JSON file.

    It is an error if the keyword does not exist or is not
    convertable to an integer.

    Parameters
    ----------
    infile : str
        The name of a file containing JSON, with the 'masterid'
        field.

    Returns
    -------
    mid : int
        The masterid value.

    """

    jcts = read_json(infile)
    try:
        mid = jcts['masterid']
    except KeyError:
        raise IOError("masterid keyword missing in {}".format(infile))

    try:
        return int(mid)
    except ValueError:
        raise IOError("unexpected masterid={} in {}".format(mid,
                                                            infile))


def read_mhulls_json(datadir, userdir, ensemble, revision,
                     validate=True):
    """Read in the master hull information (JSON).

    This reads in all the master-hull information stored in
    JSON in the datadir and userdir directories. It can also
    enforce that each master is listed as one of "accept",
    "delete", or "manual". There is no check that manual
    hulls have a polygon (a later check will do this).

    Parameters
    ----------
    datadir : str
        The location of the data directory for the ensemble.
        This directory contains the mhull files and original
        JSON files.
    userdir : str
        The location containing the user's choices for this
        ensemble (JSON files created by the QA server).
    ensemble : str
        The ensemble name.
    revision : int
        The revision number
    validate : bool, optional
        If set (the default) then the useraction for each hull is
        checked to make sure that it is in one of the accepted
        states ("accept", "delete", or "manual").

    Returns
    -------
    mhulls : dict
        The keys are the master id values. The values are the JSON
        contents as a dict. It is possible for this to be empty,
        either because there are no masters left in this ensemble
        (all have been deleted) or because there has been no
        user decision for this revision.
    """

    hulls = {}
    filename = make_hull_name_json(ensemble, None, revision)
    pat1 = os.path.join(datadir, ensemble, filename)
    for infile in glob.glob(pat1):

        # To use the chs_utils read logic, we end up reading this
        # file twice - the first time to get the master id value
        # since this is needed by the chs_utils version. An
        # alternative would be to deconstruct the file name to
        # extract the version number. Neither option is ideal,
        # and I don't want to change the chs_utils version at this
        # time.
        #
        mid = get_masterid_json(infile)
        assert mid not in hulls, mid

        cts = read_ensemble_hull_json(datadir, userdir,
                                      ensemble, mid, revision)
        if cts is None:
            raise IOError("No master-hull data from {}".format(infile))

        # What is the user decision for this master?
        # Note that there is a possibility that useraction=''
        #
        if validate:
            decision = get_user_setting(cts, 'useraction')
            if decision not in ['accept', 'delete', 'manual']:
                raise IOError("Master hull {} has ".format(mid) +
                              "decision='{}'".format(decision))

        hulls[mid] = cts

    return hulls


def read_component_json(datadir, userdir, ensemble,
                        stack, cpt, revision):
    """Read in the component-level information (JSON).

    Read in the user decision for this stack-level component. There
    is an error if a component has a master id of 0 (i.e. it was
    marked to be moved but no move has happened) or if it is
    marked as deleted (-1) and has any other masterid in it.

    Parameters
    ----------
    datadir : str
        The location of the data directory for the ensemble.
        This directory contains the mhull files and original
        JSON files.
    userdir : str
        The location containing the user's choices for this
        ensemble (JSON files created by the QA server).
    ensemble : str
        The ensemble name.
    stack : str
        The stack name
    cpt : int
        The component number.
    revision : int
        The revision number

    Returns
    -------
    mhulls : dict
        The component info, with an added field 'match_type'
        which is set to 'ambiguous', 'unambiguous', or 'deleted'
        depending on the number of master_id values
        (-1 is deleted, a single positive value is unambiguous,
        multiple positive values is ambiguous).

    Notes
    -----
    The error message may report the wrong file (e.g. user rather than
    data) depending on where the value was read from. I don't think
    it's worth fixing this.
    """

    filename = make_component_name_json(ensemble, stack, cpt,
                                        revision)
    infile = os.path.join(datadir, ensemble, filename)
    jcts = read_json(infile)
    if jcts is None:
        raise IOError("Unable to read: {}".format(infile))

    # Set up for user information.
    #
    midkey = 'master_id'
    userkeys = [midkey, 'include_in_centroid', 'lastmodified']
    for key in userkeys:
        if not setup_user_setting(jcts, key, stringval=False):
            raise IOError("Unable to set up {}".format(key))

    infile = os.path.join(userdir, ensemble, filename)
    if os.path.exists(infile):
        ucts = read_json(infile)
        for key in userkeys:
            if key in ucts:
                jcts[key]['user'] = ucts[key]

    # Validate the master id setting
    #   - must be at least 1 value
    #   - if multiple then all > 0
    #   - no 0 values
    #
    # [CHECK E]
    mids = get_user_setting(jcts, midkey)
    if len(mids) == 0:
        raise IOError("No master ids in {}".format(infile))
    minid = min(mids)
    if minid < -1:
        raise IOError("Master id={} in {}".format(minid, infile))
    if 0 in mids:
        raise IOError("{} has master_id=0".format(infile))
    if len(mids) > 1 and any([mid < 1 for mid in mids]):
        raise IOError("Ambiguous match but mid < 1 " +
                      "in {}".format(infile))

    # [CHECK F] -> not really a check, more a definition
    # [CHECK G] -> diffo
    if -1 in mids:
        match_type = 'deleted'
    elif len(mids) > 1:
        match_type = 'ambiguous'
    else:
        match_type = 'unambiguous'

    # I am not sure why I put this in here; I am deleting the check
    # since it is not valid, but this may be due to poor information
    # flow (i.e. allowing some data through in one of the files has
    # lead to the retention of unwanted data?)
    #
    # assert 'match_type' not in jcts, str(jcts)
    jcts['match_type'] = match_type
    return jcts


def read_poly_from_json(userdir, ensemble, mid, revision):
    """Read in the user-defined polygon for this master hull.

    It is an error for the polygon not to exist or to have
    less than 3 vertexes (and there must be only one).

    Parameters
    ----------
    userdir : str
        The location containing the user's choices for this
        ensemble (JSON files created by the QA server).
    ensemble : str
        The ensemble name.
    revision : int
        The revision number

    Returns
    -------
    poly : ndarray
        The RA and Dec of the polygon vertices. The polygon is
        closed, and has shape (2, nvertex).

    """

    filename = make_poly_name_json(ensemble, mid, revision)
    infile = os.path.join(userdir, filename)

    cts = read_json(infile)
    for key in ['ensemble', 'masterid', 'lastmodified', 'polygons',
                'revision']:
        if key not in cts:
            raise IOError("polygon {} missing ".format(infile) +
                          "key={}".format(key))

    polys = cts['polygons']
    npolys = len(polys)
    if npolys == 0:
        raise IOError("{} contains no polygons".format(infile))
    elif npolys > 1:
        raise IOError("{} contains {} polygons".format(infile,
                                                       npolys))

    poly = polys[0]
    for key in ['ra', 'dec']:
        if key not in poly:
            raise IOError("{} polygon missing ".format(infile) +
                          "{} key".format(key))

    ra = poly['ra']
    dec = poly['dec']
    nra = len(ra)
    ndec = len(dec)
    if nra != ndec:
        raise IOError("{} polygon ra/dec mismatch".format(infile))

    # [CHECK K]
    if nra < 3:
        raise IOError("{} polygon < 3 vertexes".format(infile))

    # [CHECK L]
    if (ra[0] != ra[-1]) or (dec[0] != dec[-1]):
        ra.append(ra[0])
        dec.append(dec[0])

    return np.vstack((ra, dec))


def read_poly_from_mhull(hulllist, masterid):
    """Read in the polygon from the master hull file.

    It is an error for the polygon not to exist or to have
    less than 3 vertexes.

    Parameters
    ----------
    hulllist : dict
        The second argument of chs_utils.read_mhull
    masterid : int
        The master id

    Returns
    -------
    poly, stack : ndarray, str
        The RA and Dec of the polygon vertices. The polygon is
        closed, and has shape (2, nvertex). The stack is the
        "base" stack for the polygon.

    """

    try:
        hull = hulllist[masterid]
    except KeyError:
        raise IOError("Unable to get masterid={}".format(masterid) +
                      " from mhull file")

    eqpos = hull['eqpos']
    if np.any(~np.isfinite(eqpos)):
        raise IOError("hull polygon is not finite! " +
                      "masterid={}".format(masterid))

    # [CHECK K]
    n = eqpos.shape[1]
    if n < 3:
        raise IOError("masterid={} ".format(masterid) +
                      "polygon < 3 vertexes")

    # [CHECK L]
    ra = eqpos[0]
    dec = eqpos[1]
    if (ra[0] != ra[-1]) or (dec[0] != dec[-1]):
        # It should be closed
        print("WARNING: mhull polygon not closed; " +
              "masterid={}".format(masterid))

        # Using np.resize automatically copies over the first
        # element into the new space, which is what we want. The
        # copies are because the arrays are slices, so need to
        # be converted into their own data.
        #
        ra = ra.copy()
        dec = dec.copy()
        npts = ra.shape + 1
        ra = np.resize(ra, npts)
        dec = np.resize(dec, npts)
        eqpos = np.vstack((ra, dec))

    return eqpos, hull['base_stk']
