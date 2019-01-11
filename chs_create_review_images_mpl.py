#!/usr/bin/env python

"""
Usage:

  ./chs_create_review_images_mpl.py chsfile outdir
      --mrgsrc3dir <dirname>
      --stkevt3dir <dirname>
      --stkfov3dir <dirname>

Aim:

Create the ensemble overview image and the per-master-hull image in
outdir using matplotlib. Also writes out JSON files for the ensemble
and each master hull.

The chsfile is the FITS table created by the script
chs_create_initial_masters.py, or by any down-stream modifications.

The following files will be written to outdir

  field.<ensemble>.v<revision>.png
  hull.<ensemble>.<master_id>.p<page_number>.v<revision>.<scale>png

  field.<ensemble>.v<revision>.json
  hull.<ensemble>.<master_id>.v<revision>.json
  cpt.<ensemble>.<stack>.<cpt>.v<revision>.json

where revision is the CHSVER keyword in chsfile, and revision,
master_id, and page_number are written out as three character,
zero-padded strings. scale is the scaling used in the image display
and can be one of log10, sqrt, none.

"""

import json
import os

from collections import defaultdict

import pycrates

import chs_status
import chs_utils as utils
import chs_review_plots_mpl as plots

# do this after import plots above, to set the backend
from matplotlib import pyplot as plt


help_str = "Create the review PNG files for CHS in this ensemble."


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


def create_review_products(chsfile, outdir,
                           mrgsrc3dir, stkevt3dir,
                           stkfov3dir, xmdat3dir):
    """Create the review products.

    Parameters
    ----------
    chsfile : str
        The FITS file containing the master hull data.
    outdir : str
        The output directory, which may be created by the routine.
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

    chsdir = os.path.dirname(chsfile)

    hullmatch, hulllist, metadata = utils.read_master_hulls(chsfile,
                                                            mrgsrc3dir)

    ensemble = metadata['ensemble']
    ensemblemap = metadata['ensemblemap']
    revision = metadata['revision']

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
        if chs_status.is_qa(src['status']):
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

    if not os.path.isdir(outdir):
        os.mkdir(outdir)

    # Context image: FOV + all the hulls
    #
    hulls = [hullmap[s] for s in stacks]
    plots.draw_ensemble_outline(ensemble, hulllist, hulls, qas,
                                fov3files)

    filename = 'field.{}.v{:03d}.png'.format(ensemble, revision)
    outfile = os.path.join(outdir, filename)

    plt.savefig(outfile)
    print("Created: {}".format(outfile))

    # Create the ensemble JSON file
    #
    revstr = "{:03d}".format(revision)
    ensdata = {'name': ensemble,
               'revision': revstr,
               'nmasters': len(mids),
               'nstacks': len(stacks),
               'ncpts': ncpts,
               'stackmap': ensemblemap,
               'status': chs_status.TODO,
               'lastmodified': '',  # could add date string here
               'usernotes': ''
               }

    filename = utils.make_field_name_json(ensemble, revision)
    outfile = os.path.join(outdir, filename)
    open(outfile, 'w').write(json.dumps(ensdata))
    print("Created: {}".format(outfile))

    # Create the per-component JSON files
    #
    for stkhulls in hullmatch.values():
        for stkhull in stkhulls:

            # Need to convert from NumPy booleans to Python ones
            # otherwise the serialization to JSON fails.
            #
            stkdata = {'lastmodified': '',
                       'stack': stkhull['stack'],
                       'component': stkhull['component'],
                       'key': stkhull['key'],
                       'ensemble': ensemble,
                       'revision': revstr,
                       # Note: need an array for master ids (in case
                       #       of ambiguous links)
                       'master_id': [stkhull['master_id']],
                       'likelihood': stkhull['likelihood'],
                       'eband': stkhull['eband'],
                       'mrg3rev': stkhull['mrg3rev'],
                       'mancode':
                       bool(stkhull['mancode']),
                       'stksvdqa':
                       bool(stkhull['stksvdqa']),
                       'include_in_centroid':
                       bool(stkhull['include_in_centroid'])}

            filename = utils.make_component_name_json(ensemble,
                                                      stkhull['stack'],
                                                      stkhull['component'],
                                                      revision)
            outfile = os.path.join(outdir, filename)
            open(outfile, 'w').write(json.dumps(stkdata))
            print("Created: {}".format(outfile))

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
            pinfo = plots.draw_hulls_and_images(src,
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

        if chs_status.is_qa(src['status']):
            action = 'manual'
        else:
            action = ''

        ensdata = {'ensemble': ensemble,
                   'masterid': mid,
                   'revision': revstr,
                   'ncpts': len(hullmatch[mid]),
                   'npages': pinfo['npages'],
                   'useraction': action,
                   'usernotes': ''
                   }
        outfile = utils.save_master(outdir, ensdata)
        print("Created: {}".format(outfile))


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

    create_review_products(args.chsfile, args.outdir,
                           mrgsrc3dir=args.mrgsrc3dir,
                           stkevt3dir=args.stkevt3dir,
                           stkfov3dir=args.stkfov3dir,
                           xmdat3dir=args.xmdat3dir)
