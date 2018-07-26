#!/usr/bin/env python

"""

Usage:

   ./chs_finalize_ensemble.py datadir userdir ensemble

Aim:

Check that the given ensemble has been through QA and is valid.
If so, create a finalized mhull file.

Checks include:

   [A] the ensemble has been marked as completed (field JSON)
   [B] the per-hull and per-component files agree with this (JSON)

   [C] all stack components from CHSVER=1 are included
   [D] no stack components are added to those in CHSVER=1
       (this is done implicitly by using the CHSVER=1 values as the
       set that is looked at)
   [E] stack level components are either deleted (master id=-1)
       or have master id > 0
   [F] if a stack level component is unambiguously matched it is
       only associated with one master id
   [G] if a stack level component is ambiguously matched it is
       associated with multiple masters

   [H] each master has at least one stack-level component
   [I] master ids are >= 1
   [J] each master has status=okay
   [K] the number of points in the polgon is > 2 and is correct
   [L] the master polygon is closed
   [M] the master polygon is not self-intersecting

   [N] master polygons are convex hulls
   [O] no master polygons overlap
   [P] each polygon has a base stack that is a member of the ensemble

   [Q] stack-level polygons only overlap the master polygons they
       are associated with (whether unambiguous or ambiguous)

   [R] each master has at least one include_in_centroid match
       (that isn't from a deleted component)

   Warnings:

   [WA] if a stack-level component is not deleted if it's
        likelihood is < 350

A master can consist of ambiguous-only matches (or match) in
this scheme. I have asked Joe if this is a problem and it isn't.

Conversions: TODO

On input, all deleted components are expected to have Master_Id=-1,
but on output we want
   a) unique negative values
   b) an offset term added (need to check Roger's requirement)

Master_Ids to be renumbered to 1 to n from whatever (positive) values
they are and then to add an offset.

"""

import glob
import os
import sys

from collections import defaultdict

import numpy as np

import pycrates

import chs_utils as utils


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

    jcts = utils.read_json(infile)
    try:
        mid = jcts['masterid']
    except KeyError:
        raise IOError("masterid keyword missing in {}".format(infile))

    try:
        return int(mid)
    except ValueError:
        raise IOError("unexpected masterid={} in {}".format(mid,
                                                            infile))


def read_mhulls_json(datadir, userdir, ensemble, revision):
    """Read in the master hull information (JSON).

    This reads in all the master-hull information stored in
    JSON in the datadir and userdir directories. It also
    enforces that each master is listed as one of "accept",
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

    Returns
    -------
    mhulls : dict
        The keys are the master id values. The values are the JSON
        contents as a dict.
    """

    hulls = {}
    filename = utils.make_hull_name_json(ensemble,
                                         None,
                                         revision)
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

        cts = utils.read_ensemble_hull_json(datadir, userdir,
                                            ensemble, mid, revision)
        if cts is None:
            raise IOError("No master-hull data from {}".format(infile))

        # What is the user decision for this master?
        decision = utils.get_user_decision(cts, 'useraction')
        if decision not in ['accept', 'delete', 'manual']:
            raise IOError("Master hull {} has ".format(mid) +
                          "decision={}".format(decision))

        hulls[mid] = cts

    return hulls


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

    filename = utils.make_poly_name_json(ensemble, mid, revision)
    infile = os.path.join(userdir, filename)

    cts = utils.read_json(infile)
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
        raise IOError("masterid={} ".format(materid) +
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

    filename = utils.make_component_name_json(ensemble, stack, cpt,
                                              revision)
    infile = os.path.join(datadir, ensemble, filename)
    jcts = utils.read_json(infile)
    if jcts is None:
        raise IOError("Unable to read: {}".format(infile))

    # Set up for user information.
    #
    midkey = 'master_id'
    userkeys = [midkey, 'include_in_centroid', 'lastmodified']
    for key in userkeys:
        if not utils.setup_user_setting(jcts, key, stringval=False):
            raise IOError("Unable to set up {}".format(key))

    infile = os.path.join(userdir, ensemble, filename)
    if os.path.exists(infile):
        ucts = utils.read_json(infile)
        for key in userkeys:
            if key in ucts:
                jcts[key]['user'] = ucts[key]

    # Validate the master id setting
    #   - must be at least 1 value
    #   - if multiple then all > 0
    #   - no 0 values
    #
    # [CHECK E]
    mids = utils.read_user_setting(jcts, midkey)
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

    assert 'match_type' not in jcts, str(jcts)
    jcts['match_type'] = match_type
    return jcts


def create_mhull(outfile, ensemble, revision, hullmd,
                 cpts, mhulls, polys,
                 compzero=0, creator=None):
    """Create the final master hull file.

    Parameters
    ----------
    outfile : str
        The name of the file to create. This must not exist.
    ensemble : str
        The ensemble name (e.g. 'ens0000100_001').
    revision : int
        The revision number of the finalized mhull file.
    hullmd : dict
        The mhull metadata read from the last revision (from
        chs_utils.read_master_hulls).
    cpts : list
        The component information.
    mhulls : dict
        The master hulls (the keys are the master id).
    polys : dict
        The master hull outlines (the keys are the master id)
    compzero : int, optional
        The component offset to apply to the component values
        when writing out the HULLMATCH block. It is expected to be
        0.
    creator : str or None
        The value for the CREATOR keyword, if set.
    """

    header = {'compzero': compzero,
              'stacks': hullmd['stacks'],
              'svdqafile': hullmd['svdqafile'],
              'centroidfile': hullmd['centroidfile']}

    hullmatch = {}
    hullmatch['master_id'] = mid
    hullmatch['nhulls'] = nhulls
    hullmatch['stackid'] = stks
    hullmatch['component'] = cpts
    hullmatch['match_type'] = mtypes
    hullmatch['area'] = areas
    hullmatch['eband'] = ebands
    hullmatch['likelihood'] = lhoods
    hullmatch['man_code'] = mancodes
    hullmatch['mrg3rev'] = revnums
    hullmatch['include_in_centroid'] = incl_centroids
    hullmatch['stksvdqa'] = stksvdqas

    hulllist = {}
    hulllist['master_id'] = mid
    hulllist['status'] = status
    hulllist['base_stk'] = base_stack
    hulllist['manmatch'] = man_match
    hulllist['manreg'] = man_reg
    hulllist['nvertex'] = nvertex
    hulllist['nstkhull'] = nstkhull
    hulllist['eqpos'] = eqpos

    utils.create_mhull_file(ensemble, revision, outfile,
                            hullmatch, hulllist, header,
                            creator=creator)

    # Check we can read in both blocks, but no real validation
    # beyond that.
    #
    chk = pycrates.CrateDataset(outfile, mode='r')
    assert chk.get_ncrates() == 3
    assert chk.get_current_crate() == 2

    def valid(idx, name, nrows):
        bl = chk.get_crate(idx)
        assert bl.name == name
        assert bl.get_nrows() == nrows

    # TODO: is len(cpts) actually correct when we have ambiguous
    # matches?
    valid(1, 'PRIMARY', 0)
    valid(2, 'HULLMATCH', len(cpts))
    valid(3, 'HULLLIST', len(mhulls))
    chk = None


def indicate_finalised(datadir, ensemble, revision):
    """Add the terminal file for this ensemble.

    Add the marker file to let the system know this ensemble has
    been finished. This has two uses:

      a) easy to spot by other software (rather than having to
         work out if an mhull file is finalised)
      b) it ensures that we can recover from a failure during
         finalization, which may leave a mhull file on disk
         that has not been validated.

    Parameters
    ----------
    datadir : str
        The location of the data directory - that is the mhull
        and original JSON files. It is also the location where
        the final mhull file is written.
    ensemble : str
        The ensemble name (e.g. 'ens0000100_001').
    revision : int
        The revision number of the finalized mhull file (this
        file must exist in datadir).

    """

    checkfile = utils.make_hull_name(ensemble, revision)
    infile = os.path.join(datadir, ensemble, checkfile)

    # Just check that the file exists (ie no validation).
    #
    if not os.path.isfile(infile):
        raise IOError("Unable to find {}".format(infile))

    outname = utils.make_terminal_name(ensemble)
    outfile = os.path.join(datadir, ensemble, outname)
    if os.path.isfile(outfile):
        raise IOError("Terminal file already exists: {}".format(outfile))

    utils.touch_file(outfile)
    if not os.path.isfile(outfile):
        # should not be needed
        raise IOError("Unable to create: {}".format(outfile))


def finalize(datadir, userdir, ensemble,
             creator=None,
             mrgsrc3dir='/data/L3/chs_master_match/input/mrgsrc3'):
    """Create the final mhull file for this ensemble.

    Parameters
    ----------
    datadir : str
        The location of the data directory - that is the mhull
        and original JSON files. It is also the location where
        the final mhull file is written.
    userdir : str
        The location from which the user has run the QA server.
    ensemble : str
        The ensemble name (e.g. 'ens0000100_001').
    creator : str or None
        If set, used for the value of the CREATOR keyword in the
        output file.
    mrgsrc3dir : str, optional
        The directory containing the mrgsrc3 files.

    Notes
    -----
    The CHSVER=1 version of the mhull file is used to define
    all stack-level components that are in this ensemble.

    """

    # Convert to the relevant ensemble-level directory
    datadir_ = os.path.join(datadir, ensemble)
    userdir_ = os.path.join(userdir, ensemble)

    for dirname in [datadir_, userdir_]:
        if not os.path.isdir(dirname):
            raise IOError("No directory called: {}".format(dirname))

    # What does the mhull say?
    mhulls = utils.find_mhulls(datadir_, ensemble)
    latest = mhulls[-1]
    revision = latest[0]
    hullname = latest[1]
    print("Ensemble {}  version {:03d}".format(ensemble, revision))

    # Read in the input master-hull data for this revision and
    # note down if there has been a change in MRGSRC3 revision.
    #
    vals = utils.read_master_hulls(hullname, mrgsrc3dir)
    hullmatch, hulllist, hullmd = vals

    # Also read in the revision=1 values so we know what the
    # total set of components is. This is done separately
    # even if the final revision is revision=1.
    #
    # It is not clear if this is a sensible step, since it
    # requires carrying along "deleted" components from previous
    # revisions. However, we need this information in the output
    # file for Joe, so whether we carry it along or add in at
    # this point is an open question.
    #
    assert mhulls[0][0] == 1, mhulls[0][0]
    vals1 = utils.read_master_hulls(mhulls[0][1], mrgsrc3dir)
    expected_components = []
    for cpts in vals1[0].values():
        # hullmatch is a dict (key = masterid) where the
        # values is a list of component info.
        #
        for cpt in cpts:
            expected_components.append((cpt['stack'],
                                        cpt['component']))

    # what is the ensemble status?
    status = utils.read_ensemble_status(datadir, userdir, ensemble,
                                        revision)

    # [CHECK A]
    if status != 'done':
        sys.stderr.write("ERROR: ensemble is marked " +
                         "status={}\n".format(status))
        sys.exit(1)

    # What are the component-level decisions?
    #
    # This ensures that all components from CHSVER=1 have
    # information and that there are no "master_id=0"
    # cases, or ones marked as both deleted and not.
    #
    # TODO: what is the best way to index this information?
    #
    # [CHECK C]
    # [CHECK D] => as mentioned this is implicit rather than explicit
    #
    cpts = []
    for stack, cpt in expected_components:
        cptinfo = read_component_json(datadir, userdir, ensemble,
                                      stack, cpt, revision)
        cpts.append(cptinfo)

    # Choice of masterid vs master_id, so use a variable
    #
    midkey = 'master_id'

    # read in the JSON for each master hull and ensure that it is
    # "done" (accepted, deleted, or manually modified).
    #
    # Those that are labelled as manually-modified are checked to
    # ensure we have a polygon for them. If not the polygon is read
    # from the input mhull file.
    #
    # For each master check that if it is valid then it has
    # components; if it is deleted then it has no components.
    #
    # [CHECK B]
    mhulls_json = read_mhulls_json(datadir, userdir, ensemble,
                                   revision)
    polys = {}
    for mid, cts in mhulls_json.items():
        decision = utils.get_user_decision(cts, 'useraction')

        # polygon information
        if decision == 'manual':
            poly = read_poly_from_json(userdir, ensemble,
                                       mid, revision)
            basestk = None

        elif decision == 'accept':
            poly, basestk = read_poly_from_mhull(hulllist, mid)

        elif decision == 'delete':
            poly = None
            basestk = None

        else:
            raise RuntimeError("Unexpected " +
                               "decision={}".format(decision))

        polys[mid] = {'eqpos': poly, 'base_stk': basestk}

        # Are there any associated components?
        has_cpt = False
        for cpt in cpts:
            mids = utils.get_user_decision(cpt, midkey)
            if mid in mids:
                has_cpt = True
                break

        if decision == 'delete':
            if has_cpt:
                raise IOError("Master id {} is ".format(mid) +
                              "deleted but it has components")
        elif not has_cpt:
            raise IOError("Master id {} has ".format(mid) +
                          "no components")

    # Ensure that the components do not refer to an unknown master.
    #
    # [CHECK H] -> kind of; this check doesn't really fit
    for cpt in cpts:
        for mid in utils.get_user_decision(cpt, midkey):
            if mid == -1:
                continue

            if mid not in mhulls_json:
                raise IOError("cpt references unknown master:\n" +
                              "{}".format(cpt))

    # Ensure that each master has at least one 'include_in_centroid'.
    #
    # Ignore any "deleted" component that has include_in_centroid
    # set to True as this shouldn't be acted on by anything.
    #
    # Note: include_in_centroid is really a per-stack setting than a
    # per-component one. So ensure that all components in a stack
    # match (at least for a given hull).
    #
    # Note that ambiguous components can have include_in_centroid
    # since it is really the stack not the component.
    #
    # This is done after the above checks (rather than being part
    # of it) as it may simplify a check slightly.
    #
    # [CHECK R]
    for mid, cts in mhulls_json.items():
        # do not care about deleted masters
        if utils.get_user_decision(cts, 'useraction') == 'delete':
            continue

        stacks = defaultdict(set)
        for cpt in cpts:
            if mid not in utils.get_user_decision(cpt, midkey):
                continue

            stack = cpt['stack']
            flag = cpt['include_in_centroid']
            stacks[stack].add(flag)

        if len(stacks) == 0:
            raise IOError("Master id {} ".format(mid) +
                          "has no include_in_centroid")

        for stack, flags in stacks.items():
            if len(flags) != 1:
                raise IOError("stack {} is in/out ".format(stack) +
                              "of master id {}".format(mid))

    """
TODO: look at polygons
  - have a base stack
  - non self-intersecting
  - do not overlap
    """

    # Ensure we can write to the directory:
    # https://stackoverflow.com/a/2113511
    #
    if not os.access(datadir, os.W_OK | os.X_OK):
        raise IOError("unable to write to {}".format(datadir))

    # QUS:
    #   what are the conditions that we do not need to create a
    #   new file? Actually, we need to do so since the mhull
    #   file we initially create is, by design, not finished
    #   (as it requires the status field to be changed).
    #
    finished_revision = revision + 1
    outname = utils.make_hull_name(ensemble, finished_revision)
    outfile = os.path.join(datadir, ensemble, outname)
    if os.path.exists(outfile):
        raise IOError("Output file exists: {}".format(outfile))

    create_mhull(outfile, ensemble, finished_revision,
                 hullmd,
                 cpts, mhulls_json, polys,
                 creator=creator)

    # have decided do not want to do a validation step after
    # reading in the file, as this repeats the checks we have
    # already done
    #
    indicate_finalised(datadir, ensemble, finished_revision)


def usage(progName):
    sys.stderr.write("Usage: {} datadir userdir " +
                     "ensemble\n".format(progName))
    sys.stderr.write("\n  datadir is the directory containing the ")
    sys.stderr.write("mhull and image files.\n")
    sys.stderr.write("  userdir is the directory from which the ")
    sys.stderr.write("QA server is run.\n")
    sys.stderr.write("  ensemble is the ensemble to review\n")
    sys.stderr.write("\n")
    sys.exit(1)


if __name__ == "__main__":

    nargs = len(sys.argv)
    if nargs != 4:
        usage(sys.argv[0])

    datadir = os.path.normpath(os.path.abspath(sys.argv[1]))
    userdir = os.path.normpath(os.path.abspath(sys.argv[2]))
    ensemble = sys.argv[3]

    finalize(datadir, userdir, ensemble,
             creator=sys.argv[0])
