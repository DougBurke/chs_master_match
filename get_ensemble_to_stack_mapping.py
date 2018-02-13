#!/usr/bin/env python

"""
Usage:
  ./get_ensemble_to_stack_mapping.py

Aim:

Query the sqlite database and create a mapping from ensemble to stack,
which is written in text format to stdout. The output format is
four columns: ensemble, nstacks, ctr, and stack. The columns - which
are space-separated - are the ensmeble name, the number of stacks in
the ensemble, the position within the current ensemble (starting at
1), and the stack name.

The ordering of the stacks within the ensemble should not be taken to
have any specific meaning (at present it matches the order retrieved
from the catops database).


"""

import sqlite3


DATABASE = '/data/L3/ops/db/catops.sqlite'


def query_db(dbase):
    """Return the ensemble/stack mapping.

    Parameters
    ----------
    dbase : str
        The name of the database to query. It is expected to contain
        the ensembles table with the rows ensemble_name and cohorts.

    Returns
    -------
    mapping : dict
        The keys are the ensemble value and the values are the
        stack contents, as a list.

    Notes
    -----
    For now there is no error checking on the db queries; there are
    some sanity checks on the returned data.

    """

    conn = sqlite3.connect(dbase)
    c = conn.cursor()

    qry = 'SELECT ensemble_name,cohorts FROM ensembles'
    out = {}
    for ens, stklist in c.execute(qry):
        if ens in out:
            raise IOError("Multiple occurences of ensemble {}".format(ens))

        stks = stklist.split(",")
        if len(stks) == 0:
            raise IOError("ensemble {} is empty".format(ens))

        def check(s):
            return s.startswith('acisf') or s.startswith('hrcf')

        if not(all([check(s) for s in stks])):
            raise IOError("Unexpected element in ensemble {}:\n{}".format(ens, stklist))

        out[ens] = stks

    return out


def doit():

    mapping = query_db(dbase=DATABASE)

    print("# ensemble nstacks ctr stack")

    for ens in sorted(mapping.keys()):

        stks = mapping[ens]
        nstk = len(stks)
        hdr = "{} {}".format(ens, nstk)

        for i, stk in enumerate(stks):
            print("{} {} {}".format(hdr, i + 1, stk))


if __name__ == "__main__":

    import sys
    if len(sys.argv) != 1:
        sys.stderr.write("Usage: {}\n".format(sys.argv[0]))
        sys.exit(1)

    doit()
