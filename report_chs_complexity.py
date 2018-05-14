#!/usr/bin/env python

"""

Usage:

  ./report_chs_complexity.py stkfile hulldir

Aim:

Report on the "complexity" of the convex-hull sources, where
complexity is the number of stacks and obsids that the hull covers.

For the moment this is ALL stacks in the ensemble, which is an upper
estimate on the complexity.

The directory is that created by the review process - i.e. where the
master-hulls are stored.

The stkfile parameter is a list of the stacks along with the number of
observations in this stack (using the columns stack and nobs).

"""

import glob
import os

from collections import defaultdict

import numpy as np
import pycrates


def read_stkfile(infile):
    """Read in the stack information.

    Parameters
    ----------
    infile : str
        A file readable by the DM with columns 'stack' and 'nobs',
        which give the number of obis in the stack.

    Returns
    -------
    ans : dict
        The keys are the stack name and the value the number of
        obis in the stack.

    """

    cr = pycrates.read_file(infile + "[cols stack, nobs]")
    out = {}
    for stack, nobs in zip(cr.stack.values,
                           cr.nobs.values.astype(np.int)):

        # try and remove the numpy string
        stack = str(stack)

        assert stack not in out, stack
        out[stack] = nobs

    return out


def read_ensemble(infile):
    """Read in the ensemble file.

    Parameters
    ----------
    infile : str
        The master-hull file for this ensemble. It must have
        HULLMATCH and HULLLIST blocks.

    Notes
    -----
    At present only HULLLIST is looked for and used.

    """

    ds = pycrates.CrateDataset(infile, mode='r')

    # We are primarily interested in the HULLLIST block, but may
    # need to look at HULLMATCH for QA cases.
    #
    cr = ds.get_crate('HULLLIST')
    nhulls = cr.get_nrows()
    if nhulls == 0:
        raise IOError("No rows in HULLLIST block of {}".format(infile))

    ver = cr.get_key_value('CHSVER')
    ensemble = cr.get_key_value('ENSEMBLE')

    stacks = []
    for i in range(cr.get_key_value('STKIDNUM')):
        stacks.append(cr.get_key_value('STKID{:03d}'.format(i)))

    if stacks == []:
        raise IOError("No stacks stored in {}".format(infile))

    coords = []
    for nv, eqpos in zip(cr.NVERTEX.values,
                         cr.EQPOS.values):
        if nv == 0:
            cs = None
        else:
            cs = eqpos[:, :nv].copy()

        coords.append(cs)

    return {'infile': infile,
            'version': ver,
            'ensemble': ensemble,
            'stacks': stacks,
            'hulls': coords}


def find_ensembles(indir):
    """Find the master-hull files.

    May report multiple files for the same ensemble, if there are
    multiple versions.
    """

    pat = os.path.join(indir, 'ens*_001', 'master_hulls.ens*fits')
    return glob.glob(pat)


def read_latest(indir):

    store = {}
    for infile in find_ensembles(indir):
        out = read_ensemble(infile)
        ens = out['ensemble']
        ver = out['version']

        if ens not in store:
            store[ens] = out
            continue

        oldver = store[ens]['version']
        if ver > oldver:
            store[ens] = out

    return store


def process_latest(stkfile, indir):

    stks = read_stkfile(stkfile)
    ensembles = read_latest(indir)

    print("# Hull directory: {}".format(indir))
    print("# Stack file: {}".format(stkfile))
    print("# Found {} ensembles with master hulls".format(len(ensembles)))

    # What is the total number of hulls?
    ntot = 0
    for ensemble in ensembles.values():
        ntot += len(ensemble['hulls'])

    print("# Total number of master hulls: {}".format(ntot))

    # Record the number of stacks and number of observations each
    # hull covers.
    #
    # At present this over-estimates the complexity since there is
    # no check if the stack covers the hull (or the obsid in a
    # stack).
    #
    nstacks = defaultdict(int)
    nobs = defaultdict(int)

    for ensemble in ensembles.values():

        nhulls = len(ensemble['hulls'])
        ns = len(ensemble['stacks'])

        nobis = sum([stks[s] for s in ensemble['stacks']])

        nstacks[ns] += nhulls
        nobs[nobis] += nhulls

    print("#")
    print("# Nstacks   Nhulls     %hulls")
    for key in sorted(list(nstacks.keys())):
        print("{:3d}  {:2d}  {:6.2f}".format(key, nstacks[key],
                                             100.0 * nstacks[key] / ntot))

    print("#")
    print("# Nobis   Nhulls     %hulls")
    for key in sorted(list(nobs.keys())):
        print("{:3d}  {:2d}  {:6.2f}".format(key, nobs[key],
                                             100.0 * nobs[key] / ntot))


if __name__ == "__main__":

    import sys
    if len(sys.argv) != 3:
        sys.stderr.write("Usage: {} stkfile hulldir\n".format(sys.argv[0]))
        sys.exit(1)

    process_latest(sys.argv[1], sys.argv[2])
