#!/usr/bin/env python

"""

Usage:

   ./chs_complete_ensemble.py datadir userdir ensemble
       --mrgsrc3dir /path/to/mrgsrc3files/
Aim:

Indicate that this ensemble has been through QA and can be
finalized - i.e. the next step would be to call

   chs_finalize_ensemble.py datadir userdir ensemble

This does not have all the checks chs_finalize_ensemble.py has,
instead it focusses on managing the deletion step (i.e. ensuring
that the data matches up).

At present it only supports 'accept' and 'reject'; 'manual' masters
cause the script to error out.

"""

import json
import os
import sys
import time

from collections import defaultdict

import numpy as np

import chs_utils as utils


def create_mhull(outfile, ensemble, revision, hullmd,
                 cpts, mhulls, polys,
                 mstzero=0,
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

    """

    assert mstzero >= 0

    header = {'compzero': compzero,
              'mstzero': mstzero,
              'stacks': hullmd['stacks'],
              'svdqafile': hullmd['svdqafile'],
              'centroidfile': hullmd['centroidfile']}

    hullmatch = defaultdict(list)
    for cpt in cpts:

        # A component can have multiple rows in this table
        mids = utils.get_user_setting(cpt, 'master_id')

        incl_in_centroid = utils.get_user_setting(cpt,
                                                  'include_in_centroid')

        if len(mids) == 1:
            assert cpt['match_type'] in ['unambiguous', 'deleted'], \
                cpt
        else:
            assert cpt['match_type'] == 'ambiguous', cpt

        for mid in mids:

            assert mid != 0, mids
            if mid < 0:
                nhulls = -1
            else:
                nhulls = mhulls[mid]['ncpts']

            hullmatch['master_id'].append(mid)

            hullmatch['nhulls'].append(nhulls)

            hullmatch['stackid'].append(cpt['stack'])
            hullmatch['component'].append(cpt['component'])
            hullmatch['match_type'].append(cpt['match_type'])

            # Do not really care about the area any longer (and am not
            # carrying the value around), so just use a terminal value.
            #
            # hullmatch['area'].append(cpt['area'])
            hullmatch['area'].append(np.nan)
            hullmatch['eband'].append(cpt['eband'])
            hullmatch['likelihood'].append(cpt['likelihood'])

            hullmatch['man_code'].append(cpt['man_code'])

            hullmatch['mrg3rev'].append(cpt['mrg3rev'])
            hullmatch['include_in_centroid'].append(incl_in_centroid)
            hullmatch['stksvdqa'].append(cpt['stksvdqa'])

    # Be explicit in the ordering here
    #
    hulllist = defaultdict(list)
    for mid in sorted(list(mhulls.keys())):

        mhull = mhulls[mid]
        status = utils.get_user_setting(mhull, 'useraction')
        if status == 'delete':
            continue

        poly = polys[mid]
        eqpos = poly['eqpos']
        nvertex = eqpos.shape[1]

        # This will have to be updated, but at present, when only
        # supporting deletion or acceptance of the proposed hull,
        # neither flag is set.
        #
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


def indicate_completed(datadir, ensemble, revision,
                       creator=None):
    """Add the field JSON file saying it is done. Maybe.

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
    creator : str, optional
        The name of the file (used in the usernotes field).

    """

    checkfile = utils.make_mhull_name(ensemble, revision)
    infile = os.path.join(datadir, ensemble, checkfile)

    # Just check that the file exists (ie no validation).
    #
    if not os.path.isfile(infile):
        raise IOError("Unable to find {}".format(infile))

    outname = utils.make_field_name_json(ensemble,
                                         revision=revision)
    outfile = os.path.join(datadir, ensemble, outname)
    if os.path.isfile(outfile):
        raise IOError("Status=done file already exists: {}".format(outfile))

    out = {'status': 'done',
           'lastmodified': time.asctime(),
           'revision': "{:03d}".format(int(revision))}

    if creator is None:
        out['usernotes'] = 'auto completed'
    else:
        out['usernotes'] = 'completed by {}'.format(creator)

    with open(outfile, 'w') as fh:
        fh.write(json.dumps(out))


def complete(datadir, userdir, ensemble,
             creator=None,
             mrgsrc3dir='/data/L3/chs_master_match/input/mrgsrc3'):
    """Create the finished mhull file for this ensemble.

    This assumes that the hull has either been deleted or
    the CHSVER=1 version accepted. It does not YET support
    modifying the hull or re-generating the master hull from
    a subset of the original components.

    Parameters
    ----------
    datadir : str
        The location of the data directory - that is the mhull
        and original JSON files. It is also the location where
        the completed mhull file is written.
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
    # master_ids_base has the Master_Id as the key and
    # contains a set of (stack, component) values, representing
    # the original set of stack-level hulls that contributed
    # to the master.
    #
    assert mhulls[0][0] == 1, mhulls[0][0]
    vals1 = utils.read_master_hulls(mhulls[0][1], mrgsrc3dir)
    expected_components = []
    master_ids_base = defaultdict(set)

    original_component_data = defaultdict(list)

    # hullmatch (first argument of read_master_hulls) is a
    # dict (with key = masterid) where the values are a
    # list of component info.
    #
    for mid, cpts in vals1[0].items():

        assert mid not in master_ids_base

        for cpt in cpts:
            key = (cpt['stack'], cpt['component'])
            master_ids_base[mid].add(key)
            expected_components.append(key)

            # components can be associated with multiple masters,
            # so store each entry in a list (this is wasteful since
            # all the other fields are the same).
            #
            original_component_data[key].append(cpt)

    # what is the ensemble status?
    status = utils.read_ensemble_status(datadir, userdir,
                                        ensemble, revision)

    # This script converts review to done: it may get changed
    # to require a "done" setting, which would be set by the
    # UI, but currently we do not have this capability.
    #
    if status != 'review':
        sys.stderr.write("ERROR: ensemble is marked " +
                         "status={}\n".format(status))
        sys.exit(1)

    # Choice of masterid vs master_id, so use a variable
    #
    midkey = 'master_id'

    # Read in the component-level decisions (based on the ver=001
    # list of components).
    #
    # There is limited checking of how consistent everything is.
    #
    # master_ids_new is master_ids_base but for this version.
    # This will be modified when reviewing the master-hull
    # decisions (if the hull is marked as being deleted).
    #
    cpts = []
    master_ids_new = defaultdict(set)
    for key in expected_components:
        stack, cpt = key
        cptinfo = utils.read_component_json(datadir, userdir,
                                            ensemble, stack, cpt,
                                            revision)

        # The mancode value in the JSON file is a boolean, but
        # we want the actual integer value from the input. This
        # is messy.
        #
        odata = original_component_data[key]
        mancodes = set([o['mancode'] for o in odata])
        man_codes = set([o['man_code'] for o in odata])

        assert len(mancodes) == 1, mancodes
        assert len(man_codes) == 1, man_codes

        mancode = mancodes.pop()
        man_code = man_codes.pop()

        assert mancode == cptinfo['mancode']

        assert 'man_code' not in cptinfo
        cptinfo['man_code'] = man_code

        cpts.append(cptinfo)

        for mid in utils.get_user_setting(cptinfo, midkey):
            assert mid != 0
            if mid < 0:
                continue

            master_ids_new[mid].add(key)

    # read in the JSON for each master hull and ensure that it is
    # "done" (accepted, deleted, or manually modified).
    #
    # Those that are labelled as manually-modified are checked to
    # ensure we have a polygon for them. If not the polygon is read
    # from the input mhull file.
    #
    # For each master check that if it is valid then it has
    # components; if it is deleted then ensure any components
    # it has are removed from this master.
    #
    # [CHECK B]
    #
    # NOTE: any manually-modified hulls cause the script to error
    #       out at this time (as this is not supported in the first
    #       released version).
    #
    mhulls_json = utils.read_mhulls_json(datadir, userdir, ensemble,
                                         revision)
    polys = {}
    for mid, cts in mhulls_json.items():

        assert mid in master_ids_new, mid

        decision = utils.get_user_setting(cts, 'useraction')

        # polygon information
        if decision == 'manual':

            # This is a temporary check, just for the first pass of sources
            emsg = "Master_Id={} ".format(mid) + \
                "has useraction={} ".format(decision) + \
                "which is not supported in this version\n"
            sys.stderr.write(emsg)
            sys.exit(1)

            poly = utils.read_poly_from_json(userdir, ensemble,
                                             mid, revision)
            basestk = None

        elif decision == 'accept':
            poly, basestk = utils.read_poly_from_mhull(hulllist, mid)

        elif decision == 'delete':
            poly = None
            basestk = None

            # component deletion handled below
        else:
            raise RuntimeError("Unexpected " +
                               "decision={}".format(decision))

        polys[mid] = {'eqpos': poly, 'base_stk': basestk}

        if decision == 'delete':

            del master_ids_new[mid]

            # Could use master_ids_new[mid] to find the components,
            # but would then have to loop through cpts anyway as do
            # not have a better indexing system at this time.
            #
            for cpt in cpts:
                mids = utils.get_user_setting(cpt, midkey)
                if mid not in mids:
                    continue

                # Restrict the supported cases for now
                if len(mids) != 1:
                    sys.stderr.write("ERROR: component deletion from ambiguous match not supported for master={}\n".format(mid))
                    sys.exit(1)

                mids.remove(mid)
                nmids = len(mids)
                if nmids == 0:
                    mids = [-1]
                    cpt['match_type'] = 'deleted'

                elif nmids == 1:
                    assert mids[0] > 0, mids

                    # component is now unambiguously matches
                    cpt['match_type'] = 'unambiguous'

                cpt[midkey]['user'] = mids

    # For the first round, ensure that the master ids from
    # CHSVER=1 are either deleted or all accepted (that is,
    # the components are all deleted or all accepted).
    #
    for mid1, cpts1 in master_ids_base.items():

        # What is the decision for this master?
        cts = mhulls_json[mid1]
        decision = utils.get_user_setting(cts, 'useraction')

        if decision == 'delete':
            assert mid1 not in master_ids_new, mid1
            continue

        assert cpts1 == master_ids_new[mid1], mid1

    if not utils.have_directory_write_access(datadir):
        raise IOError("unable to write to {}".format(datadir))

    finished_revision = revision + 1
    outname = utils.make_mhull_name(ensemble, finished_revision)
    outfile = os.path.join(datadir, ensemble, outname)
    if os.path.exists(outfile):
        raise IOError("Output file exists: {}".format(outfile))

    create_mhull(outfile, ensemble, finished_revision,
                 hullmd,
                 cpts, mhulls_json, polys,
                 creator=creator)

    indicate_completed(datadir, ensemble, finished_revision,
                       creator=creator)


help_str = "Mark an ensemble as done (ready for chs_finalize_ensemble)."


if __name__ == "__main__":

    import argparse

    parser = argparse.ArgumentParser(description=help_str,
                                     prog=sys.argv[0])

    parser.add_argument("datadir",
                        help="The directory containing the mhull and image files")
    parser.add_argument("userdir",
                        help="The directory from which the QA server is run")
    parser.add_argument("ensemble", type=str,
                        help="The ensemble to complete")

    parser.add_argument("--mrgsrc3dir",
                        default="/data/L3/chs_master_match/input/mrgsrc3",
                        help="The mrgsrc3 directory: default %(default)s")

    args = parser.parse_args(sys.argv[1:])

    datadir = os.path.normpath(os.path.abspath(args.datadir))
    userdir = os.path.normpath(os.path.abspath(args.userdir))

    complete(datadir, userdir, args.ensemble,
             mrgsrc3dir=args.mrgsrc3dir,
             creator=sys.argv[0])
