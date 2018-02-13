#!/usr/bin/env python

"""
Usage:

  ./get_ensembles_with_hulls.py ensemblelist stacklist

Aim:

From the ensemble list, which is assumed to have columns 'ensemble'
and 'stack' and so is used to map between stack and ensemble
identifier, and the stacklist, which contains a list of stacks
(one per line) which contain stack-level hulls (i.e. MEXTSRC block
has at least one row with STATUS=0), output a list of those
ensembles which contain hulls.

The output is to the screen, one ensemble per line.

"""

import pycrates
import stk


def doit(ensemblename, stackname):
    """Find ensembles with hulls.

    The ensemble names are displayed to the screen, in ascending
    order, one per line (so usable as a stack).

    Parameters
    ----------
    ensemblename : str
        The name of a file which contains columns ensemble and
        stack, and so gives the mapping from stack to ensemble.
    stackname : str
        The name of an ASCII file which has one stack per line,
        and contains stacks with hulls.

    """

    cr = pycrates.read_file(ensemblename + "[cols ensemble,stack]")
    stackmap = {}
    for ensname, stkname in zip(cr.ensemble.values, cr.stack.values):
        assert stkname not in stackmap, stkname
        stackmap[stkname] = ensname

    cr = None

    ensembles = set([])
    for stkname in stk.build("@" + stackname):
        try:
            ensembles.add(stackmap[stkname])
        except KeyError:
            raise IOError("Unrecognized stack '{}'".format(stkname))

    if len(ensembles) == 0:
        raise IOError("No ensembles found!")

    for ensname in sorted(list(ensembles)):
        print(ensname)


if __name__ == "__main__":

    import sys
    if len(sys.argv) != 3:
        sys.stderr.write("Usage: {} ".format(sys.argv[0]) +
                         "ensemblelist stacklist\n")
        sys.exit(1)

    doit(sys.argv[1], sys.argv[2])

