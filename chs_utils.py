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

    return sorted(list(zip(vers, matches)),
                  key=lambda x: x[0])


def make_field_name_json(ensemble, revision):
    """What is the JSON file storing the field data?

    Parameters
    ----------
    ensemble : str
        The ensemble name.
    revision : int
        The revision number.

    Returns
    -------
    filename : str
        The name of the JSON file (with no path).
    """

    # assert isinstance(revision, six.string_types)
    revision = int(revision)

    return 'field.{}.v{:03d}.json'.format(ensemble, revision)


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
        store = {'master_id': mid,
                 'stack': stackid,
                 'component': component - compzero,
                 'key': key,
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
    case the 'proposed' setting is used.

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

    if proposed is None or proposed == '':
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
    status : {'todo', 'review', 'done'}

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

    if ans is None:
        warn("status=None for ensemble " +
             "{} version {}".format(ensemble, revision))
        return "unknown"

    errlog("no status for ensemble " +
           "{} version {}".format(ensemble, revision))
    return "unknown"


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

    Parameters
    ----------
    userdir : str
        The location for the user-stored data.
    data : dict
        The JSON dictionary containing the elements to write out.
    """

    ensemble = data['name']
    version = data['revision']  # this is in string form, 0 padded

    outdir = os.path.join(userdir, ensemble)
    if os.path.exists(outdir):
        if not os.path.isdir(outdir):
            raise IOError("Exists but not a directory! {}".format(outdir))
    else:
        os.mkdir(outdir)

    outname = make_field_name_json(ensemble, version)
    outfile = os.path.join(outdir, outname)

    store = {"name": ensemble,
             "lastmodified": time.asctime(),
             "usernotes": data['usernotes'],
             "status": data['status'],
             "revision": version}
    with open(outfile, 'w') as fh:
        fh.write(json.dumps(store))


def save_master(userdir, data):
    """Save the master-level details.

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

    outname = make_hull_name_json(ensemble, masterid, version)
    outfile = os.path.join(outdir, outname)

    store = {"ensemble": ensemble,
             "masterid": masterid,
             "lastmodified": time.asctime(),
             "usernotes": data['usernotes'],
             "useraction": data['useraction'],
             "revision": version}
    with open(outfile, 'w') as fh:
        fh.write(json.dumps(store))


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

    Parameters
    ----------
    userdir : str
        The location for the user-stored data.
    data : dict
        The JSON dictionary containing the elements to write out.
    """

    ensemble = data['ensemble']
    version = data['revision']  # this is in string form, 0 padded
    stack = data['stack']
    component = data['component']

    outdir = os.path.join(userdir, ensemble)
    if os.path.exists(outdir):
        if not os.path.isdir(outdir):
            raise IOError("Exists but not a directory! {}".format(outdir))
    else:
        os.mkdir(outdir)

    outname = make_component_name_json(ensemble,
                                       stack,
                                       component,
                                       version)
    outfile = os.path.join(outdir, outname)

    # For now we do not copy over the other fields (that should be
    # read only) that are in the "original" version of this file
    # (such as likelihood, band, mrg3rev).
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
    with open(outfile, 'w') as fh:
        fh.write(json.dumps(store))


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
