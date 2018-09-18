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

"""

import os
import sys

from collections import defaultdict

import pycrates

import chs_utils as utils


def create_mhull(outfile, ensemble, revision, hullmd,
                 cpts, mhulls, polys,
                 mstzero=9000,
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
    mstzero : int, optional
        The offset to apply to the Master_Id values for both the
        HULLMATCH and HULLLIST blocks. It must be 0 or greater.
        This value is written out as the MSTZERO keyword in the
        header.
    compzero : int, optional
        The component offset to apply to the component values
        when writing out the HULLMATCH block. It is expected to be
        0. This value is written out as the COMPZERO keyword
        in the header.
    creator : str or None
        The value for the CREATOR keyword, if set.

    Notes
    -----
    Deleted components are re-numbered to have unique master ids:
    -1, -2, ...
    """

    """DBG:
    print("***"); print(hullmd)
    print("***"); print(cpts)
    print("***"); print(mhulls)
    print("***"); print(polys)
    print(polys[1]['eqpos'].shape[1])
    """

    assert mstzero >= 0

    header = {'compzero': compzero,
              'mstzero': mstzero,
              'stacks': hullmd['stacks'],
              'svdqafile': hullmd['svdqafile'],
              'centroidfile': hullmd['centroidfile']}

    # track the deleted components
    ndel = 0

    hullmatch = defaultdict(list)
    for cpt in cpts:

        # A component can have multiple rows in this table
        mids = utils.get_user_setting(cpt, 'master_id')

        incl_in_centroid = utils.get_user_setting(cpt,
                                                  'include_in_centroid')

        # TODO: check that this hasn't already been done,
        #       because it should have
        #
        #       Should this support delete?
        #
        if len(mids) == 1:
            assert cpt['match_type'] == 'unambiguous', cpt
        else:
            assert cpt['match_type'] == 'ambiguous', cpt

        for mid in mids:

            assert mid != 0, mids
            if mid < 0:
                ndel += 1

                mid = -1 * ndel
                nhulls = -1

            else:
                nhulls = mhulls[mid]['ncpts']

            hullmatch['master_id'].append(mid)

            hullmatch['nhulls'].append(nhulls)

            hullmatch['stackid'].append(cpt['stack'])
            hullmatch['component'].append(cpt['component'])
            hullmatch['match_type'].append(cpt['match_type'])

            # TODO: why is area info missing in input?
            hullmatch['area'].append(-1.0)
            hullmatch['eband'].append(cpt['eband'])
            hullmatch['likelihood'].append(cpt['likelihood'])

            # TODO: has man_code been corrected?
            hullmatch['man_code'].append(cpt['mancode'])

            hullmatch['mrg3rev'].append(cpt['mrg3rev'])
            hullmatch['include_in_centroid'].append(incl_in_centroid)
            hullmatch['stksvdqa'].append(cpt['stksvdqa'])

    # Be explicit in the ordering here
    #
    hulllist = defaultdict(list)
    for mid in sorted(list(mhulls.keys())):

        mhull = mhulls[mid]
        poly = polys[mid]
        eqpos = poly['eqpos']
        nvertex = eqpos.shape[1]

        status = utils.get_user_setting(mhull, 'useraction')

        # TODO: fix this
        man_match = False
        man_reg = False

        hulllist['master_id'].append(mid)
        hulllist['status'].append(status)
        hulllist['base_stk'].append(poly['base_stk'])
        hulllist['manmatch'].append(man_match)
        hulllist['manreg'].append(man_reg)
        hulllist['nvertex'].append(nvertex)
        hulllist['nstkhull'].append(mhull['ncpts'])
        hulllist['eqpos'].append(eqpos)

    utils.create_mhull_file(ensemble, revision, outfile,
                            hullmatch, hulllist, header,
                            creator=creator)

    # Check we can read in both blocks, but no real validation
    # beyond that.
    #
    chk = pycrates.CrateDataset(outfile, mode='r')
    assert chk.get_ncrates() == 3
    assert chk.get_current_crate() == 2

    # Note: hack for nrows=None to mean empty table
    def valid(idx, name, nrows=None):
        bl = chk.get_crate(idx)
        assert bl.name == name
        if nrows is None:
            assert bl.get_shape() == (0,)
        else:
            assert bl.get_nrows() == nrows

    # TODO: is len(cpts) actually correct when we have ambiguous
    # matches?
    valid(1, 'PRIMARY')
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

    checkfile = utils.make_mhull_name(ensemble, revision)
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
    status = utils.read_ensemble_status(datadir, userdir,
                                        ensemble, revision)

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
        cptinfo = utils.read_component_json(datadir, userdir,
                                            ensemble, stack, cpt,
                                            revision)
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
    mhulls_json = utils.read_mhulls_json(datadir, userdir, ensemble,
                                         revision)
    polys = {}
    for mid, cts in mhulls_json.items():
        decision = utils.get_user_setting(cts, 'useraction')

        # polygon information
        if decision == 'manual':
            poly = utils.read_poly_from_json(userdir, ensemble,
                                             mid, revision)
            basestk = None

        elif decision == 'accept':
            poly, basestk = utils.read_poly_from_mhull(hulllist, mid)

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
            mids = utils.get_user_setting(cpt, midkey)

            print("DBG: mid = {}  mids = {}".format(mid, mids))

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
        for mid in utils.get_user_setting(cpt, midkey):
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
        if utils.get_user_setting(cts, 'useraction') == 'delete':
            continue

        stacks = defaultdict(set)
        for cpt in cpts:
            if mid not in utils.get_user_setting(cpt, midkey):
                continue

            stack = cpt['stack']
            flag = utils.get_user_setting(cpt, 'include_in_centroid')
            # print("DBG: {} -> {}".format(cpt['include_in_centroid'], flag))
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

    if not utils.have_directory_write_access(datadir):
        raise IOError("unable to write to {}".format(datadir))

    # QUS:
    #   what are the conditions that we do not need to create a
    #   new file? Actually, we need to do so since the mhull
    #   file we initially create is, by design, not finished
    #   (as it requires the status field to be changed).
    #
    finished_revision = revision + 1
    outname = utils.make_mhull_name(ensemble, finished_revision)
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
    sys.stderr.write("Usage: {} datadir userdir ".format(progName) +
                     "ensemble\n")
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
