#!/usr/bin/env python

"""
Usage:

  ./chs_create_review_images.py chsfile outdir
      --mrgsrc3dir <dirname>
      --stkevt3dir <dirname>
      --stkfov3dir <dirname>

Aim:

Create the ensemble overview image and the per-master-hull image in
outdir.

The chsfile is the FITS table created by the script
chs_create_initial_masters.py, or by any down-stream modifications.

The following files will be written to outdir

  field.<ensemble>.v<revision>.png
  hull.<ensemble>.<master_id>.p<page_number>.v<revision>.<scale>png

where revision is the CHSVER keyword in chsfile, and revision,
master_id, and page_number are written out as three character,
zero-padded strings. scale is the scaling used in the image display
and can be one of log10, sqrt, none.

"""

import os

import pycrates
import pychips.all as pychips

import chs_utils as utils
import chs_review_plots as plots

help_str = "Create the review PNG files for CHS in this ensemble."


def read_qa_hulls(qadir, revision, master_id):
    """Read in the QA hulls for the given master hull.

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
        (closed, only finite values).

    """

    infile = os.path.join(qadir,
                          'qa.{:03d}.v{:03d}.fits'.format(master_id,
                                                          revision))
    cr = pycrates.read_file(infile)
    if cr.get_nrows() < 2:
        raise IOError("Expected at least 2 rows")

    out = []
    for cpt, npts, eqpos in zip(cr.COMPONENT.values,
                                cr.NVERTEX.values,
                                cr.EQPOS.values):
        eqpos = utils.validate_polygon(eqpos[:, :npts], report=True)
        out.append({'component': cpt,
                    'master_id': master_id,
                    'eqpos': eqpos.copy()})

    return out


def create_review_products(chsfile, outdir,
                           mrgsrc3dir, stkevt3dir, stkfov3dir,
                           showhack=False):
    """Create the review products.

    Parameters
    ----------
    chsfile : str
        The FITS file containing the master hull data.
    outdir : str
        The output directory, which may be created by the routine.
    mrgsrc3dir, stkevt3dir, stkfov3dir : str
        The directory names containing the mrgsrc3, evt3, and fov3
        files for the stacks. The names must match
        <stack>*<type>.fits[.gz] and there can only be one per stack
        per type.
    showhack : bool, optional
        Does the window have to be displayed before it can be printed?
        This seems to depend on the system.

    Notes
    -----
    Should the "context" plot - i.e. all hulls from an ensemble - also
    include the current master hulls? I worry that it will make it
    harder to see the details, but let's you see at a glance if there
    are potential issues.
    """

    chsdir = os.path.dirname(chsfile)

    srcmatch, srclist, metadata = utils.read_master_hulls(chsfile)

    ensemble = metadata['ensemble']
    revision = metadata['revision']

    # What stacks do we care about (those with hulls)
    stacks = set([])

    for stkhulls in srcmatch.values():
        for stkhull in stkhulls:
            stacks.add(stkhull['stack'])

    assert len(stacks) > 0

    stacks = list(stacks)

    # It might be more useful to sort the stacks by "number of fov
    # files" say, or some other criteria.
    #
    stacks = sorted(stacks)

    mrgsrc3files = [utils.find_mrgsrc3(s, mrgsrc3dir)
                    for s in stacks]
    fov3files = [utils.find_stkfov3(s, stkfov3dir)
                 for s in stacks]

    if not os.path.isdir(outdir):
        os.mkdir(outdir)

    # Context image: FOV + all the hulls
    #
    hullmap = {s: utils.read_hulls_from_mrgsrc3(m)
               for s, m in zip(stacks, mrgsrc3files)}

    hulls = [hullmap[s] for s in stacks]
    plots.draw_ensemble_outline(ensemble, hulls, fov3files)

    filename = 'field.{}.v{:03d}.png'.format(ensemble, revision)
    outfile = os.path.join(outdir, filename)

    if showhack:
        pychips.set_window(['display', True])

    pychips.print_window(outfile, ['clobber', True])
    print("Created: {}".format(outfile))

    # Per-master hulls
    #
    for mid in sorted(srclist.keys()):

        src = srclist[mid]
        if src['status'] == 'qa':
            qahulls = read_qa_hulls(chsdir, revision,
                                    src['master_id'])
        else:
            qahulls = None

        # a lot of repeated work to support different scalings,
        # but not worth the complexity of avoiding this
        for evtscale in ['log10', 'sqrt', 'none']:
            plots.draw_hulls_and_images(src,
                                        srcmatch[mid],
                                        hullmap,
                                        stkevt3dir,
                                        outdir,
                                        ensemble, revision,
                                        evtscale=evtscale,
                                        qahulls=qahulls,
                                        showhack=showhack)


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
                        default="/data/L3/jmiller/archive_mrgsrc3",
                        help="The mrgsrc3 directory: default %(default)s")
    parser.add_argument("--stkevt3dir", type=str,
                        default="/data/dburke2/L3/rel2.0/master_match/convex_hulls/stkevt3",
                        help="The stkevt3 directory: default %(default)s")
    parser.add_argument("--stkfov3dir", type=str,
                        default="/data/dburke2/L3/rel2.0/master_match/convex_hulls/stkfov3",
                        help="The stkfov3 directory: default %(default)s")
    parser.add_argument("--showhack", action='store_true',
                        help="Display ChIPS window before printing?")

    args = parser.parse_args(sys.argv[1:])

    create_review_products(args.chsfile, args.outdir,
                           mrgsrc3dir=args.mrgsrc3dir,
                           stkevt3dir=args.stkevt3dir,
                           stkfov3dir=args.stkfov3dir,
                           showhack=args.showhack)
