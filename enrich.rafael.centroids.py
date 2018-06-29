#!/usr/bin/env python

"""

Usage:

  ./enrich.rafael.centroids.py indir infile

Aim:

Given Rafael's list of masters and the choice of HRC or ACIS for the
centroid calculation, output a list of "stack component include"
lines to the screen.

indir should point to the directory containing the master-hull data;
which should be something like
    /data/L3/chs_master_match/test/2018-03-09a/

infile is Rafael's input, which is expected to be the file
rafael.centroids.2018-06-14.dat

For those master hulls which we do not have info - at present this
is Cen A since the image data was not available for Rafael to review -
we mark all the components as being exclude which will make the user
choose.

"""

from collections import defaultdict
import os

import pycrates

import chs_utils as utils


# do not want to rely on the now-standard Python Enum or its backport
# so using https://www.pythoncentral.io/how-to-implement-an-enum-in-python/
#
def enum(*args):
    enums = dict(zip(args, range(len(args))))
    return type('Enum', (), enums)


# ACIS only, HRC only, mixed, or user has to decide
#
# Not 100% clear to me at present on the intended meaning of the
# manual tag in Rafael's email. Actually, looking at the text again
# I believe the MANUAL tag is not relevant for this process.
#
Option = enum('UNKNOWN', 'ACIS', 'HRC', 'MIXED')


def parse_rafael(infile):
    """Parse Rafael's file."""

    store = defaultdict(list)
    ctr = 0

    with open(infile, 'r') as fh:
        for l in fh.readlines():
            if l.startswith('#'):
                continue

            # toks = l.strip().split(maxsplit=4)  Python 3-ism
            toks = l.strip().split(None, 4)
            ens = toks[1]
            mid = int(toks[2])
            rhs = toks[4]

            # assume only two tokens; let it error out if not
            nacis, nhrc = [int(t) for t in toks[3].split(',')]

            # the "decision" column is normally a single token but can
            # be "hrc_only, man"
            #
            toks = rhs.split()
            opt = toks[0]

            if opt.endswith(','):
                if toks[1] != "man":
                    raise ValueError(l.strip)

                # For now we ignore the manual flag
                # choice = option.MANUAL
                opt = opt[:-1]

            if opt == 'ign_hrc':
                choice = Option.ACIS
            elif opt == 'hrc_only':
                choice = Option.HRC
            elif opt == 'inc_hrc':
                choice = Option.MIXED
            elif opt == 'Cas':
                # This is a hack
                choice = Option.UNKNOWN
            else:
                raise IOError("Unexpected option={}".format(opt))

            store[ens].append({'masterid': mid, 'nacis': nacis,
                               'nhrc': nhrc, 'choice': choice})
            ctr += 1

    assert ctr == 18, ctr
    return store


def read_ensemble(indir, ensemble):
    """Read master-hull file."""

    pat = os.path.join(indir, ensemble,
                       utils.make_mhull_name(ensemble))
    match = utils.find_single_match(pat)

    infile = "{}[HULLMATCH][cols Master_Id,STACKID,COMPONENT]".format(match)
    cr = pycrates.read_file(infile)
    compzero = cr.get_key_value('COMPZERO')
    if compzero is None:
        raise IOError("No COMPZERO in {}".format(infile))

    store = defaultdict(list)
    for mid, stack, cpt in zip(cr.Master_Id.values,
                               cr.STACKID.values,
                               cr.COMPONENT.values):
        store[mid].append((stack, cpt - compzero))

    out = {}
    for mid, cpts in store.items():
        nacis = len([s for s in cpts if s[0].startswith('acis')])
        nhrc = len([s for s in cpts if s[0].startswith('hrc')])
        assert nacis + nhrc == len(cpts)

        out[mid] = {'ensemble': ensemble,
                    'components': cpts, 'nacis': nacis, 'nhrc': nhrc}

    return out


def print_choice(hull, choice):
    """Is the given hull included or excluded from the centroid?

    Parameters
    ----------
    hull : (str, int)
        The stackid and component number.
    choice : Option enumeration

    """

    stack, cpt = hull

    is_acis = stack.startswith('acis')
    if choice == Option.ACIS:
        flag = is_acis
    elif choice == Option.HRC:
        flag = not is_acis
    elif choice == Option.MIXED:
        flag = True
    elif choice == Option.UNKNOWN:
        flag = False
    else:
        raise RuntimeError("Unexpected choice={}".format(choice))

    print("{:24s} {:2d} {}".format(stack, cpt, flag))


def find_components(mids, hulls):
    """Match Rafael's hulls to the masters and output the stacks.

    We primarily match on master id, but also do a check on nacis,nhrc
    as a secondary check. There is no attempt to look for the actual
    value master id if they do not match.
    """

    # map from masterid in hulls (Rafael's reported value) to the
    # master ids in mids. This should be an identity mapping; in fact
    # it is an error if it is not.
    #
    # The error cases are not expected to happen, so not really
    # providing enough information to track down the problem.
    #
    for hull in hulls:
        mid = hull['masterid']
        nacis = hull['nacis']
        nhrc = hull['nhrc']

        try:
            store = mids[mid]
        except KeyError:
            raise IOError("No master id match for {}".format(hull))

        if store['nacis'] != nacis or store['nhrc'] != nhrc:
            raise IOError("Nacis/Nhrc mismatch for {}".format(hull))

        for cpt in store['components']:
            print_choice(cpt, hull['choice'])


def print_header():

    print("# stack cpt use_cen")


def doit(indir, infile):
    """Process Rafael's file.

    The input format is rather hard-coded to the expected format.
    """

    raf = parse_rafael(infile)

    print_header()

    store = {}
    for ens, hulls in raf.items():

        try:
            mids = store[ens]
        except KeyError:
            store[ens] = read_ensemble(indir, ens)
            mids = store[ens]

        find_components(mids, hulls)


if __name__ == "__main__":

    import sys
    if len(sys.argv) != 3:
        sys.stderr.write("Usage: {} indir infile\n".format(sys.argv[0]))
        sys.exit(1)

    doit(sys.argv[1], sys.argv[2])
