#!/usr/bin/env python

"""

Usage:

   ./chs_complete_ensemble.py datadir userdir ensemble

Aim:

Indicate that this ensemble has been through QA and can be
finalized - i.e. the next step would be to call

   chs_finalize_ensemble.py datadir userdir ensemble

This does not have all the checks chs_finalize_ensemble.py has,
instead it focusses on managing the deletion step (i.e. ensuring
that the data matches up).

"""

import json
import os
import sys
import time

from collections import defaultdict

import chs_utils as utils


hackField = True
print("\n*** WARNING hackField hack is in operation ***\n\n")


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
        status = utils.get_user_setting(mhull, 'useraction')
        if status == 'delete':
            continue

        poly = polys[mid]
        eqpos = poly['eqpos']
        nvertex = eqpos.shape[1]

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
    if hackField:
        status = utils.read_ensemble_status(datadir,
                                            # userdir,
                                            "/pool7/dburke/test",
                                            ensemble, revision)
    else:
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

    # Read in the component-level decisions (based on the ver=001
    # list of components).
    #
    # These can be changed by the master decisions (primarily
    # delete). There is limited checking of how consistent
    # everything is.
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
    # components; if it is deleted then ensure any components
    # it has are marked as deleted.
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

        if decision == 'delete':

            for cpt in cpts:
                mids = utils.get_user_setting(cpt, midkey)
                if mid not in mids:
                    continue

                mids.remove(mid)
                if len(mids) == 0:
                    mids = [-1]
                    cpt['match_type'] = 'deleted'

                cpt[midkey]['user'] = mids

    # START HACK
    print("\n*** HACKING OUTDIR:")
    datadir = '/pool7/dburke/hulls'
    dname = os.path.join(datadir, ensemble)
    if not os.path.exists(dname):
        os.mkdir(dname)

    # END HACK

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

    complete(datadir, userdir, ensemble,
             creator=sys.argv[0])
