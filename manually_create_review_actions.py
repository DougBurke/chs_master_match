#!/usr/bin/env python

"""
Usage:

  ./manually_create_review_actions.py datadir userdir

Aim:

Runs though all the ensembles in userdir and adds in the necessary
field.ensemble.<version>.json files to indicate that the ensemble
has been through the initial version of review in the UI. This
automates the "Finish ensemble" button on the ensemble page of the
UI.

"""

import glob
import os

import chs_utils as utils
import chs_status


def find_ensembles(userdir):
    """Return all sub-directories matching ens???????_001 that
    do not have a field.ens???????_001.*.json file

    This misses out ensembles that are completed but have an outdated
    field file, but rather be safe than sorry. It is also not clever
    enough to ignore field fiels with an old version number.

    For now do not check that the ?'s are numeric/valid

    """

    pat = os.path.join(userdir, 'ens???????_001')
    out = []
    for m in glob.glob(pat):

        if not os.path.isdir(m):
            print("Warning: skipping {} as not a directory".format(m))
            continue

        # should use an os.path function
        toks = [t for t in os.path.split(m) if t != '']
        ensemble = toks[-1]

        # look for any version; ideally would ignore old versions
        # but not worth coding that at this time
        #
        fieldname = utils.make_field_name_json(ensemble)
        matches = glob.glob(os.path.join(m, fieldname))
        if len(matches) > 0:
            print("Skipping {} as has field JSON file".format(ensemble))
            continue

        out.append(ensemble)

    return sorted(out)


def compare(datadir, userdir, ensemble):
    """Add a field.ens<>.v<>.json file if it is valid to do so.

    Returns
    -------
    flag : bool
        This is true if a field file was created in userdir.

    Notes
    -----
    This assumes that there are no user-notes for the ensemble. As
    the user notes are stored in the field file, and this should
    only being run on ensembles with no field files, this assumption
    should be okay. A "dummy" user note is added to say the field
    file was auto-generated.
    """

    # The following is based on chs_finalize_ensemble.py, but greatly
    # simplified.
    #

    datadir_ = os.path.join(datadir, ensemble)
    # userdir_ = os.path.join(userdir, ensemble)

    # What does the mhull say?
    mhulls = utils.find_mhulls(datadir_, ensemble)
    latest = mhulls[-1]
    revision = latest[0]
    # hullname = latest[1]

    revstr = "{:03d}".format(revision)
    print("Ensemble {}  version {}".format(ensemble, revstr))

    # What is the ensemble status? This assumes that the decision is
    # coming from the data directory; there is currently no way to
    # enforce this constraint (the assumption here is that we have
    # checked for no field directory in userdir, so any field status
    # must be coming from the datadir).
    #
    status = utils.read_ensemble_status(datadir, userdir, ensemble,
                                        revision)

    if status != chs_status.TODO:
        print("   - skipping as status={}".format(status))
        return False

    # Check that all hulls have a decision: note this is a
    # very-limited check
    #
    mhulls = utils.read_mhulls_json(datadir, userdir, ensemble,
                                    revision, validate=False)
    for mid, mhull in mhulls.items():
        decision = utils.get_user_setting(mhull, 'useraction')
        if decision not in ['accept', 'delete', 'manual']:
            print("  - skipping as master={} ".format(mid) +
                  "has action='{}'".format(decision))
            return False

    # At this point we can create a field file
    #
    data = {'name': ensemble,
            'revision': revstr,
            'status': chs_status.REVIEW,
            'usernotes': 'The review status was auto-generated.'
            }
    utils.save_ensemble(userdir, data)
    return True


def process(datadir, userdir):
    """Look for all ensembles in userdir and add field.json files if needed."""

    ensembles = find_ensembles(userdir)
    if len(ensembles) == 0:
        print("No missing field files found in {}".format(userdir))
        return

    print("Looking at {} ensembles".format(len(ensembles)))
    ndone = 0
    for ensemble in ensembles:
        if compare(datadir, userdir, ensemble):
            ndone += 1

    print("")
    if ndone == 0:
        print("No field files were created.")
    elif ndone == 1:
        print("Added one field file.")
    else:
        print("Added {} field files.".format(ndone))


if __name__ == "__main__":

    import sys
    if len(sys.argv) != 3:
        sys.stderr.write("Usage: {} ".format(sys.argv[0]) +
                         "datadir userdir")
        sys.exit(1)

    datadir = os.path.normpath(sys.argv[1])
    userdir = os.path.normpath(sys.argv[2])
    for d in [datadir, userdir]:
        if not os.path.isdir(d):
            sys.stderr.write("Error: {} is not a directory\n".format(d))
            sys.exit(1)

    if not utils.have_directory_write_access(userdir):
        sys.stderr.write("Error: unable to write to " +
                         "{}\n".format(userdir))
        sys.exit(1)

    process(datadir, userdir)
