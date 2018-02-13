#!/usr/bin/env python

"""
Usage:

  ./get_ensembles_with_hulls.py ensemblelist stacklist

Aim:

From the ensemble list, which is assumed to have columns 'ensemble'
and 'stack' and so is used to map between stack and ensemble
identifier, and the stacklist, which contains the name of the
mrgsrc3 file (including trailing .gz but no directory specifier)
and the number of valid stack-level hulls (i.e. MEXTSRC block
has at least one row with STATUS=0), output those ensembles with
hulls, along with the number of stack-level hulls in that stack.

"""

from collections import defaultdict

import pycrates
import stk


def doit(ensemblename, stackname):
    """Find ensembles with hulls.

    The ensemble names are displayed to the screen, in ascending
    order, one per line, along with the total number of stack-level
    hulls in the ensemble (they may overlap or not).

    Parameters
    ----------
    ensemblename : str
        The name of a file which contains columns ensemble and
        stack, and so gives the mapping from stack to ensemble.
    stackname : str
        The name of an ASCII file which has two columns: the
        name of the mrgsrc3 file (no directory) and the number
        of valid hulls in that stack.
        and contains stacks with hulls.

    """

    cr = pycrates.read_file(ensemblename + "[cols ensemble,stack]")
    stackmap = {}
    for ensname, stkname in zip(cr.ensemble.values, cr.stack.values):
        assert stkname not in stackmap, stkname
        stackmap[stkname] = ensname

    cr = None

    ensembles = defaultdict(lambda: 0)
    cr = pycrates.read_file(stackname)
    for mfile, nhulls in zip(cr.get_column(0).values,
                             cr.get_column(1).values.astype(int)):
        if nhulls < 1:
            continue

        idx = mfile.find('N')
        if idx == -1:
            raise IOError("Invalid mrgsrc3 file name " +
                          "'{}'".format(mfile))

        stackname = mfile[:idx]
        try:
            ensemblename = stackmap[stackname]
        except KeyError:
            raise IOError("Unrecognized stack from " +
                          "'{}'".format(mfile))

        ensembles[ensemblename] += nhulls

    if len(ensembles) == 0:
        raise IOError("No ensembles found!")

    print("# ensemble nstackhulls")
    for ensname in sorted(list(ensembles.keys())):
        print("{} {}".format(ensname, ensembles[ensname]))


if __name__ == "__main__":

    import sys
    if len(sys.argv) != 3:
        sys.stderr.write("Usage: {} ".format(sys.argv[0]) +
                         "ensemblelist stacklist\n")
        sys.exit(1)

    doit(sys.argv[1], sys.argv[2])
