#!/usr/bin/env python

"""
Usage:

  ./chs_update_review_images.py chsfile outdir
      --mrgsrc3dir <dirname>
      --stkevt3dir <dirname>
      --stkfov3dir <dirname>

Aim:

Create the ensemble overview image and the per-master-hull image in
outdir using matplotlib using the user's updated specifications (e.g.
deleting a component or master hull).

The following files will be written to outdir

  field.<ensemble>.v<revision>.png
  hull.<ensemble>.<master_id>.p<page_number>.v<revision>.<scale>png

Unlike chs_create_review_images_mpl.py it does *NOT* create any
json files.

"""

import os

from collections import defaultdict

import pycrates

import chs_utils as utils
import chs_review_plots_mpl as plots

# do this after import plots above, to set the backend
from matplotlib import pyplot as plt


help_str = "Update the review PNG files for CHS in this ensemble."


def read_qa_hulls(qadir, revision, master_id):
    """Read in the QA hull (or hulls) for the given master hull.

    Parameters
    ----------
    qadir : str
        The directory containing the qa.<master_id>.v<revision>.fits
        files to read.
    revision : int
        The revision number to use.
    master_id : int
        The master id.

    Returns
    -------
    qahulls : list of dict
        Each entry contains the 'eqpos' keyword, which is a 2 by npts
        NumPy array with the celestial coordinates of the polygon
        (closed, only finite values). There can be one or more items.

    """

    infile = os.path.join(qadir,
                          'qa.{:03d}.v{:03d}.fits'.format(master_id,
                                                          revision))
    cr = pycrates.read_file(infile)
    if cr.get_nrows() < 1:
        raise IOError("Expected at least 1 row: {}".format(infile))

    out = []
    for cpt, npts, eqpos in zip(cr.COMPONENT.values,
                                cr.NVERTEX.values,
                                cr.EQPOS.values):
        eqpos = utils.validate_polygon(eqpos[:, :npts], report=True)
        out.append({'component': cpt,
                    'master_id': master_id,
                    'eqpos': eqpos.copy()})

    return out


def update_review_products(chsfile, outdir,
                           mrgsrc3dir, stkevt3dir,
                           stkfov3dir, xmdat3dir):
    """Update the review products.

    Parameters
    ----------
    chsfile : str
        The FITS file containing the master hull data.
    outdir : str
        The output directory, which must exist.
    mrgsrc3dir, stkevt3dir, stkfov3dir, xmdat3dir : str
        The directory names containing the mrgsrc3, evt3, fov3, and
        xmdat3 files for the stacks. The names must match
        <stack>*<type>.fits[.gz] and there can only be one per stack
        per type. The xmdat3 files are optional and are stored as
        <stack>/<stack>N000_xmdat3.fits.

    Notes
    -----
    Should the "context" plot - i.e. all hulls from an ensemble - also
    include the current master hulls? I worry that it will make it
    harder to see the details, but let's you see at a glance if there
    are potential issues.
    """

    # Since we are assumed to be running on updated files, the output
    # directory should already exist.
    #
    if not os.path.isdir(outdir):
        sys.stderr.write("ERROR: outdir does not exist ")
        sys.stderr.write("{}\n".format(outdir))
        sys.exit(1)

    chsdir = os.path.dirname(chsfile)

    hullmatch, hulllist, metadata = utils.read_master_hulls(chsfile,
                                                            mrgsrc3dir)

    ensemble = metadata['ensemble']
    ensemblemap = metadata['ensemblemap']
    revision = metadata['revision']

    filename = 'field.{}.v{:03d}.png'.format(ensemble, revision)
    ctxfile = os.path.join(outdir, filename)
    if os.path.exists(ctxfile):
        sys.stderr.write("ERROR: context image already exists ")
        sys.stderr.write("{}\n".format(ctxfile))
        sys.exit(1)

    mids = sorted(hulllist.keys())

    # What stacks do we care about (those with hulls)
    stacks = set([])

    # ncpts is the number of stack hulls in the ensemble
    ncpts = 0
    for stkhulls in hullmatch.values():
        for stkhull in stkhulls:
            stacks.add(stkhull['stack'])
            ncpts += 1

    assert len(stacks) > 0

    # It might be more useful to sort the stacks by "number of fov
    # files" say, or some other criteria.
    #
    stacks = list(stacks)
    stacks = sorted(stacks)

    # Force a check that we can find these files before any
    # processing.
    #
    fov3files = [utils.find_stkfov3(s, stkfov3dir) for s in stacks]

    # Read in all the QAs, labelled by mid.
    #
    qas = {}
    for mid in mids:
        src = hulllist[mid]
        if src['status'].startswith('qa'):
            qas[mid] = read_qa_hulls(chsdir, revision,
                                     src['master_id'])

    # Ideally I would remove hullmap, and use hullmatch directly,
    # but the format is slightly different, so recreate hullmap
    # (it used to be formed by calling read_hulls_from_mrgsrc3
    # but this is now part of read_master_hulls).
    #
    hullmap = defaultdict(list)
    for masterhull in hullmatch.values():
        for hull in masterhull:
            stack = hull['stack']
            hullmap[stack].append(hull)

    # Context image: FOV + all the hulls
    #
    hulls = [hullmap[s] for s in stacks]
    plots.draw_ensemble_outline(ensemble, hulllist, hulls, qas,
                                fov3files)

    plt.savefig(ctxfile)
    print("Created: {}".format(ctxfile))

    # Per-master hulls
    #
    for mid in mids:

        src = hulllist[mid]
        try:
            qahulls = qas[mid]
        except KeyError:
            qahulls = None

        # a lot of repeated work to support different scalings,
        # but not worth the complexity of avoiding this
        for evtscale in ['log10', 'sqrt', 'none']:
            plots.draw_hulls_and_images(src,
                                        hullmatch[mid],
                                        hullmap,
                                        stkevt3dir,
                                        xmdat3dir,
                                        outdir,
                                        ensemble,
                                        ensemblemap,
                                        revision,
                                        evtscale=evtscale,
                                        qahulls=qahulls)


if __name__ == "__main__":

    import argparse
    import sys

    parser = argparse.ArgumentParser(description=help_str,
                                     prog=sys.argv[0])

    parser.add_argument("chsfile", type=str,
                        help="The master hull data")
    parser.add_argument("outdir", type=str,
                        help="The output directory (may exist)")

    parser.add_argument("--mrgsrc3dir", type=str,
                        default="/data/L3/chs_master_match/input/mrgsrc3",
                        help="The mrgsrc3 directory: default %(default)s")
    parser.add_argument("--stkevt3dir", type=str,
                        default="/data/L3/chs_master_match/input/stkevt3",
                        help="The stkevt3 directory: default %(default)s")
    parser.add_argument("--stkfov3dir", type=str,
                        default="/data/L3/chs_master_match/input/stkfov3",
                        help="The stkfov3 directory: default %(default)s")
    parser.add_argument("--xmdat3dir", type=str,
                        default="/data/L3/chs_master_match/input/xmdat3",
                        help="The xmdat3 directory: default %(default)s")

    args = parser.parse_args(sys.argv[1:])

    update_review_products(args.chsfile, args.outdir,
                           mrgsrc3dir=args.mrgsrc3dir,
                           stkevt3dir=args.stkevt3dir,
                           stkfov3dir=args.stkfov3dir,
                           xmdat3dir=args.xmdat3dir)
