"""
Create review plots.

This requires matplotlib and CIAO.

Use
  /data/L3/chs_master_match/local/binpy27/ciao-4.9/
  /data/L3/chs_master_match/local/binpy27/ciao-4.10/

"""

from collections import defaultdict
import os

import numpy as np

import six

import pycrates
import region

import chs_utils as utils
import chs_status

from astropy.wcs import WCS
from astropy.io.fits import Header

from cycler import cycler

import matplotlib as mpl
mpl.use('Agg', warn=False)

import matplotlib.pyplot as plt
from matplotlib import colors
from matplotlib.patches import Ellipse


def tr2wcs(tr_eqpos, tr_physical=None):
    """Convert a pyTransform object to AstroPy WCS object.

    This is rather hacky, and relies on only having to worry
    about simple (e.g. WCS-TAN) transforms. It has been
    created by trial and error, rather than research, dur
    to time constraints.

    Hulls - which are tabular - just have a single transform
    but images have both sky sky and eqpos, so need hacking.
    """

    ctypes = tr_eqpos.get_parameter('CTYPE').get_value()
    crvals = tr_eqpos.get_parameter('CRVAL').get_value()
    crpixs = tr_eqpos.get_parameter('CRPIX').get_value()
    cdelts = tr_eqpos.get_parameter('CDELT').get_value()

    if tr_physical is not None:
        if tr_physical.get_parameter('ROTATION').get_value() != 0.0:
            raise ValueError("Unexpected transform")

        # convert the origin
        crpixs = tr_physical.invert([crpixs])[0]

        # convert the scale
        scale = tr_physical.get_parameter('SCALE').get_value()
        cdelts *= scale

    hdr = Header()
    hdr['CTYPE1'] = ctypes[0]
    hdr['CRVAL1'] = crvals[0]
    hdr['CRPIX1'] = crpixs[0]
    hdr['CDELT1'] = cdelts[0]
    hdr['CTYPE2'] = ctypes[1]
    hdr['CRVAL2'] = crvals[1]
    hdr['CRPIX2'] = crpixs[1]
    hdr['CDELT2'] = cdelts[1]

    return WCS(hdr)


def filter_image(vals, x, y, poly):
    """Return the list of vals which lie within the polygon.

    Parameters
    ----------
    vals, x, y : NumPy array
        The values, x coordinate and y coordinate values of
        the points to filter.
    poly : NumPy array
        2D array, (2 by npts) representing the polygon's vertices
        in the same coordinate system as x and y. The polygon
        is assumed to be closed and contain no infinite values.

    Returns
    -------
    filtered : NumPy array
        Those elements of vals that fall within the polygon.
        ACTUALLY, it now returns 0 for a pixel that was filtered out.

    """

    # would be nice to use the CIAO 4.10 region library here, so that
    # can avoid the need to create a region definition to parse.
    #
    coords = [str(v) for v in poly.T.flatten()]
    rstr = "polygon({})".format(",".join(coords))
    reg = region.regParse(rstr)

    # avoid a memory leak by looping through individual coordinate
    # pairs (I believe there is a CIAO 4.9 era memory leak if the
    # vectorized version is used).
    #
    out = []
    for (xv, yv), v in zip(np.vstack((x, y)).T, vals):
        if region.regInsideRegion(reg, xv, yv):
            out.append(v)
        else:
            out.append(0)

    return np.asarray(out)


def read_xmdat3(xmdat3dir, stack):
    """Read xmdat3 PSF data, if available.

    The xmdat3 files are assumed to be stored as
       <xmdat3dir>/<stack>/<stack>N000_xmdat3.fits

    Should I send in the ra/dec limits so that we can filter
    on these (would miss those which are centered outside the
    range but overlap)?

    If the file is missing a message is logged.

    Parameters
    ----------
    xmdat3dir : str
        The directory containing the files.
    stack : str
        The stack name

    Returns
    -------
    rval : list of regions or None
        None if there is no xmdat3 file, otherwise a list
        (which can be empty) of regions, which contain
        'ra', 'dec', 'r0', 'r1', 'angle' (ra and dec are in
        decimal degrees, r0 and r1 are in arcseconds, and angle
        is in degrees).

    """

    infile = os.path.join(xmdat3dir, stack,
                          '{}N000_xmdat3.fits'.format(stack))
    if not os.path.isfile(infile):
        utils.logonce("no xmdat3 file {}".format(infile))
        return None

    try:
        cr = pycrates.read_file(infile +
                                "[cols ra,dec,psf_r0,psf_r1,psf_ang]")
    except IOError:
        utils.logonce("unable to read XMDAT3 file {}".format(infile))
        raise

    out = []
    for ra, dec, r0, r1, ang in zip(cr.ra.values,
                                    cr.dec.values,
                                    cr.psf_r0.values,
                                    cr.psf_r1.values,
                                    cr.psf_ang.values):
        out.append({'ra': ra, 'dec': dec,
                    'r0': r0, 'r1': r1, 'angle': ang})

    return out


def draw_xmdat3(regs, axis, transform,
                edgecolor='white',
                facecolor='none'):
    """Draw on the PSF ellipses.

    This approximates the ellipse, assuming that we are on a
    tangent-plane projection and close to the tangent point.

    Parameters
    ----------
    regs : output of read_xmdat3
        If None then nothing is done.
    axis
        The axis to add the patches to.
    transform
        The WCS axis transform.
    edgecolor, facecolor : str
        Passed through to the ellipse.

    """

    if regs is None:
        return

    for r in regs:
        # Need to convert radii into diameters and then from arcsec
        # into degrees.
        #
        e = Ellipse(xy=(r['ra'], r['dec']),
                    width=r['r0'] / 1800.0, height=r['r1'] / 1800.0,
                    angle=-r['angle'],
                    transform=transform,
                    facecolor=facecolor, edgecolor=edgecolor)
        axis.add_patch(e)


# TODO: should we draw on the FOV?
#
def draw_hull(ax_trans, hull, hcolor, lthick, lstyle, opacity=1):
    "line-style handling is messy"

    if 'mancode' in hull:
        if hull['mancode'] == 0:
            lstyle = 'solid'
        else:
            lstyle = 'dotted'

    eqsrc = hull['eqpos']
    x = eqsrc[0]
    y = eqsrc[1]

    plt.plot(x, y, linestyle=lstyle, linewidth=lthick,
             alpha=opacity,
             color=hcolor, transform=ax_trans)


def label_hull(ax_trans, hull, label, color='black'):
    """Label the hull at its mid-point (well, line pointing to it)."""

    ra = hull['eqpos'][0]
    dec = hull['eqpos'][1]

    # Offset from the center to avoid a clash where the number can be
    # obscured by the hull; not sure what the best heuristic here.
    #
    ra0 = (ra.min() + ra.max()) / 2.0
    dec0 = (dec.min() + dec.max()) / 2.0

    # pick ~ 45" in both X and Y as a guess, since ~ 1' looks to
    # be the size below which confusion happens (on one test
    # case, ens0020500_001. This turned out to be a bit small
    # so increase it.
    #
    # Add in an arrow to help indicate the center of the hull.
    #
    dra = 3 * 45.0 / 3600 / np.cos(dec0 * np.pi / 180.0)
    ddec = 3 * 45.0 / 3600

    """
    I had trouble getting this to work; probably because the
    shift is small enough that the line ended up being "invisible"
    csys = ax_trans
    plt.annotate("{}".format(mid),
                 xy=(ra0, dec0),
                 xytext=(ra0 + dra, dec0 + ddec),
                 arrowprops=dict(facecolor='black',
                                 arrowstyle='-'),
                 color='black',
                 fontsize=12,
                 xycoords=csys,
                 textcoords=csys)
                 """

    plt.plot([ra0, ra0 + dra], [dec0, dec0 + ddec],
             '-', color=color,
             alpha=0.4,
             transform=ax_trans)

    plt.text(ra0 + dra, dec0 + ddec, label,
             horizontalalignment='center',
             # verticalalignment='center',
             color=color,
             fontsize=12,
             transform=ax_trans)


def find_hull_limits(stackhulls, hulldata, axscale=0.5):
    """Return limits of the hulls.

    Parameters
    ----------
    stackhulls : list
        What stack-level hulls form this master hull? Each entry
        is a dictionary with the keys 'stack' and 'component'.
    hulldata : dict
        The hull data indexed by (stack, component).
    axscale : float, optional
        The additional space around the hull, as a fraction of the
        width/height of the hull. This value refers to the delta
        added to each side (so twice this is added overall). Note
        that the space added will be more than this fraction for the
        "smaller" side of the hull, as the plot aspect ratio is
        maintained.

    Returns
    -------
    lims : dict
        The keys are: 'dmfilter', ...

    """

    ra_lims = []
    dec_lims = []
    for shull in stackhulls:
        key = shull['stack'], shull['component']
        hull = hulldata[key]
        eqsrc = hull['eqpos']
        ra = eqsrc[0]
        dec = eqsrc[1]

        ra_lims.extend([np.nanmin(ra), np.nanmax(ra)])
        dec_lims.extend([np.nanmin(dec), np.nanmax(dec)])

    # What are the limits of the hulls?
    #
    ra_lims = np.asarray(ra_lims)
    dec_lims = np.asarray(dec_lims)

    ra_min = ra_lims.min()
    ra_max = ra_lims.max()
    dec_min = dec_lims.min()
    dec_max = dec_lims.max()

    # expand limits on each side
    #
    dra = ra_max - ra_min
    ddec = dec_max - dec_min

    # Enforce a minimum padding distance (which is then
    # scaled by axscale)
    #
    minsep = 45.0 / 3600
    ddec = max(ddec, minsep)

    # divide or multiply by cos term here?
    minsep /= np.cos((dec_min + dec_max) * np.pi / (2.0 * 180.0))
    dra = max(dra, minsep)

    ra_min -= axscale * dra
    ra_max += axscale * dra

    dec_min -= axscale * ddec
    dec_max += axscale * ddec

    # Select just the region of interest to try and save memory.
    #
    sfilt = utils.make_spatial_filter_range(ra_min, ra_max,
                                            dec_min, dec_max)

    return {'dmfilter': sfilt,
            'ra': [ra_min, ra_max],
            'dec': [dec_min, dec_max]}


def label_plot(stack, cpt, bname, axes,
               hullmap, ensemblemap,
               max_lhood=500,
               lblcol='orange',
               lblsize=12):
    """Add labels to a plot of the stack.

    Parameters
    ----------
    stack : str
        The stack name.
    cpt
        The component number for the hull.
    bname : str
        The band name for the data.
    axes
        The plot axes.
    hullmap : dict
        The stack-level hull data, stored by the stack id. The values
        are those returned by read_hulls_from_mrgsrc3.
    ensemblemap : dict
        The keys are stack ids, and the values are the STKIDxxx
        value (i.e. the integer value of xxx).
    max_lhood : number, optional
        Likelihoods greater than this are displayed as "L>max_lhood"
        rather than "L=lhood"

    """

    # Convert from stack name to STKIDxxx value and add in
    # the component number.
    #
    stacklbl = "{:03d}.{:02d}".format(ensemblemap[stack], cpt)
    # bandlbl = "band={}".format(bname)
    bandlbl = bname

    xl = 0.05
    xr = 0.95
    yb = 0.05
    yt = 0.9

    plt.text(xl, yt, bandlbl,
             color=lblcol, fontsize=lblsize,
             transform=axes.transAxes)

    plt.text(xr, yb, stacklbl, horizontalalignment='right',
             color=lblcol, fontsize=lblsize,
             transform=axes.transAxes)

    # What was the likelihood?
    # Label if the stack-level hull was manually-modified;
    # ideally this would use the knowledge of whether the stack
    # went to SVDQA but for now do not pass that in (would either
    # need to be done separately or to read the mhull and not
    # mrgsrc3 file).
    #
    lhood = None
    manadj = None
    stksvdqa = None
    for hull in hullmap[stack]:
        if hull['component'] != cpt:
            continue

        assert manadj is None
        manadj = hull['mancode'] != 0

        assert lhood is None
        lhood = hull['likelihood']

        assert stksvdqa is None
        stksvdqa = hull['stksvdqa']

        # should probably break here but leave in as a sanity
        # check (that hullmap doesn't have duplicate entries)

    assert manadj is not None
    assert lhood is not None
    assert stksvdqa is not None

    if lhood > max_lhood:
        # lbl = r"$\mathcal{{L}}>{}$".format(max_lhood)
        lbl = ">{}".format(max_lhood)
    else:
        # lbl = r"$\mathcal{{L}}={:.1f}$".format(lhood)
        lbl = "{:.1f}".format(lhood)
    plt.text(xl, yb, lbl,
             color=lblcol, fontsize=lblsize,
             transform=axes.transAxes)

    lbl = None
    if manadj:
        # if manually adjusted then we don't care if it went to
        # SVD QA or not.
        # lbl = "ManAdj"
        lbl = "Man"
    elif stksvdqa:
        # lbl = "STK SVD QA"
        lbl = "SVDQA"

    if lbl is not None:
        plt.text(xr, yt, lbl, horizontalalignment='right',
                 color=lblcol, fontsize=lblsize,
                 transform=axes.transAxes)


def add_hulls_from_other_stacks(stack, stacks, hullmap, axes,
                                color='white'):
    """Draw on stack hulls from other stacks in this ensemble.

    Parameters
    ----------
    stack : str
        The stack being shown (so hulls from this stack are skipped).
    stacks : list of str
        All the stacks in the ensemble which have hulls.
    hullmap : dict
        The stack-level hull data, stored by the stack id. The values
        are those returned by read_hulls_from_mrgsrc3.
    axes
        The plot axes.
    color : str, optional
        The color to draw the hulls.
    """

    ax_trans = axes.get_transform('world')
    for ostack in stacks:
        if ostack == stack:
            continue

        for hull in hullmap[ostack]:
            draw_hull(ax_trans, hull, color, 1, 'dotted',
                      opacity=0.2)


def add_hulls(stack, cpt, hullmap, master_hull, qahulls, axes,
              master_color='gold',
              qa_color='cyan'):
    """Draw the hulls for this stack.

    draw the hulls for this stack

    draw the other hulls first, as reference
    (i.e. those from other stacks)

    NOTE: this draws on all hulls, so can be useful
          if nearby ones overlap.

    Parameters
    ----------
    stack : str
        The stack name.
    cpt
        The component number for the hull.
    hullmap : dict
        The stack-level hull data, stored by the stack id. The values
        are those returned by read_hulls_from_mrgsrc3.
    master_hull: dict
        Contains the master hull: fields are 'master_id', 'status',
        and 'eqpos'. The 'status' field should be one of:
        'todo', 'okay', 'qa[-...]'.
    qahulls : None or list of dict, optional
        This is only used if master_hull['status'] is set to 'qa[-...]'.
        Each entry represents a hull, and has the 'eqpos' field
        which contains the polygon.
    axes
        The plot axes.
    master_color : str, optional
        The color for the master hull.
    qa_color : str, optional
        The color for any QA hulls.
    """

    ax_trans = axes.get_transform('world')

    # Draw the hull for this component in orange and the
    # others in the stack as a red-ish color. Trying to
    # match masterhull.js behavior.
    #
    for hull in hullmap[stack]:
        if hull['component'] == cpt:
            hullcol = 'orange'
        else:
            hullcol = '#cc3333'

        draw_hull(ax_trans, hull, hullcol, 2, 'solid')

    # NOTE: these are drawn thinner than the stack-level hulls so
    # they do not obscure them (for cases when the two contours
    # are the same or very similar).
    #
    if chs_status.is_qa(master_hull['status']):
        for qahull in qahulls:
            draw_hull(ax_trans, qahull, qa_color, 1, 'dashed')
    else:
        draw_hull(ax_trans, master_hull, master_color, 1, 'solid')


def labelize_plot(axes, add_labels):
    """Do we add or hide plot labels?

    Parameters
    ----------
    axes
        The plot axes.
    add_labels : bool
        If True then the axes are labelled, otherwise they
        are hidden.

    """

    if add_labels:
        axes.coords['ra'].set_major_formatter('hh:mm:ss')
        axes.coords['dec'].set_major_formatter('dd:mm:ss')
    else:
        for l in ['ra', 'dec']:
            axes.coords[l].set_ticks_visible(False)
            axes.coords[l].set_ticklabel_visible(False)
            axes.coords[l].set_axislabel('')


def draw_hulls_and_images(master_hull,
                          stackhulls,
                          hullmap,
                          stkevt3dir,
                          xmdat3dir,
                          outdir,
                          ensemble,
                          ensemblemap,
                          revision,
                          qahulls=None,
                          evtscale='log10',
                          master_color='gold',
                          qa_color='cyan',
                          axscale=0.5,
                          show_other_stack_hulls=False,
                          colorbar=True):
    """Draw individidual hull+image for all members of an ensemble.

    An image is created per stack-hull; this means that if there are
    2 hulls from the same stack then there will be 2 images. This is
    a change to the original implementation, party necessitated by
    the desire to show the correct band. The display is limited to
    a "zoomed" in section around the master hull (controlled
    by the axscale). Other hulls from the stack are drawn, but will
    only be visible if they fall in this area.

    Parameters
    ----------
    master_hull: dict
        Contains the master hull: fields are 'master_id', 'status',
        and 'eqpos'. The 'status' field should be one of:
        'todo', 'okay', 'qa[-...]'.
    stackhulls : list
        What stack-level hulls form this master hull? Each entry
        is a dictionary with the keys 'stack' and 'component'.
    hullmap : dict
        The stack-level hull data, stored by the stack id. The values
        are those returned by read_hulls_from_mrgsrc3.
    stkevt3dir : str
        The location of the stkevt3 files.
    xmdat3dir : str
        The location of the xmdat3 files (per stack, optional).
    outdir : str
        The output directory, which must exist.
    ensemble : str
        The ensemble value.
    ensemblemap : dict
        The keys are stack ids, and the values are the STKIDxxx
        value (i.e. the integer value of xxx).
    revision : int
        The revision number of the file.
    qahulls : None or list of dict, optional
        This is only used if master_hull['status'] is set to 'qa[-...]'.
        Each entry represents a hull, and has the 'eqpos' field
        which contains the polygon.
    evtscale : {'none', 'log10', 'log', 'sqrt'}, optional
        The scaling to apply to the event data before displaying it.
    master_color : str, optional
        The color for the master hull.
    qa_color : str, optional
        The color for any QA hulls.
    axscale : float, optional
        The additional space around the hull, as a fraction of the
        width/height of the hull. This value refers to the delta
        added to each side (so twice this is added overall). Note
        that the space added will be more than this fraction for the
        "smaller" side of the hull, as the plot aspect ratio is
        maintained.
    show_other_stack_hulls : bool, optional
        Should the hulls from other stacks be drawn on the image?
    colorbar : bool, optional
        Should each image have a color bar?

    Notes
    -----
    For the initial revision we do not need to worry about Match_Type
    being ambiguous (this is the 'match_type' keyword of each entry
    in stackhulls). If we re-run after some master-match work then
    it may be necessary to disambiguate visually.

    """

    if not os.path.isdir(outdir):
        raise ValueError("outdir={} is not a directory".format(outdir))

    masterid = master_hull['master_id']

    # It is easier if we can map the stack hull data via the
    # key: stack, cpt. Actually, the following code is very unclear
    # about data ownership, and really could do with a clean up.
    #
    # It is also useful to know what stacks we have.
    #
    stacks = set()
    for shull in stackhulls:
        stacks.add(shull['stack'])

    stacks = sorted(list(stacks))

    # How many hulls are there in this master hull?
    nhulls = len(stackhulls)
    nstacks = len(stacks)

    # Extract all the hulls from the stacks we are interested in,
    # whether they are in related to the given master hull or not.
    #
    hulldata = {}
    for stack, hulls in six.iteritems(hullmap):
        if stack not in stacks:
            continue

        for hull in hulls:
            assert hull['stack'] == stack
            key = stack, hull['component']
            assert key not in hulldata
            hulldata[key] = hull

    # Do we have any xmdat3 files for these stacks? We want to
    # store the None to mark stacks as having no data.
    #
    xmdat3map = {}
    for stack in stacks:
        xmdat3map[stack] = read_xmdat3(xmdat3dir, stack)

    # Need limits for determining how much of the event file
    # to read in.
    #
    # Technically the QA and master hulls should be included here, but
    # they should all, by definition, be no larger than the
    # stack-level hulls.
    #
    data_lims = find_hull_limits(stackhulls, hulldata,
                                 axscale=axscale)
    sfilt = data_lims['dmfilter']

    # Find the event files for the stacks. Note that we check
    # all files and then report errors at the end (so you don't
    # have to download a file, run, download a file, run, ...
    # to get them all). Of course, if there are multiple hulls
    # in the ensemble this only handles one of them.
    #
    evtfiles = {}
    failed = []
    for stack in stacks:
        path = os.path.join(stkevt3dir,
                            '{}*_evt3.fits*'.format(stack))
        try:
            evtfiles[stack] = utils.find_single_match(path)
        except IOError as exc:
            failed.append((stack, str(exc)))

    if len(failed) > 0:
        emsg = "\n".join([f[1] for f in failed])
        raise IOError("The following event file(s) are " +
                      "missing:\n{}".format(emsg))

    # could go for a more-adaptive scheme
    if nhulls <= 9:
        nsize = np.ceil(np.sqrt(nhulls)).astype(np.int)
        nplots = nsize * nsize
        npages = 1
    else:
        nsize = 3
        nplots = nsize * nsize
        npages = np.ceil(nhulls / (nplots * 1.0)).astype(np.int)

    def save_plot(pagenum):
        """Create PNG output then destroy the window"""

        fmt = 'hull.{}.{:03d}.p{:03d}.v{:03d}.{}.png'
        outfile = fmt.format(ensemble, masterid, pagenum, revision,
                             evtscale)
        outfile = os.path.join(outdir, outfile)

        # Hide any warnings about square root of a negative number.
        #
        oldsettings = np.seterr(all='ignore')
        try:
            plt.savefig(outfile)
        finally:
            np.seterr(**oldsettings)

        print("Created: {}".format(outfile))
        plt.close()

    title = "{} {:03d}".format(ensemble, masterid) + \
        "  #stacks={} #hulls= {}".format(nstacks, nhulls)
    tcol = 'k'
    qatitle = None
    qacol = 'r'

    if chs_status.is_qa(master_hull['status']):
        qatitle = master_hull['status'].upper()

    page_idx = 0
    nplots_in_page = None
    fig = None

    lblsize = 12

    def sqrtwrapper(vmin=None, vmax=None, clip=False):
        return colors.PowerNorm(gamma=0.5, vmin=vmin, vmax=vmax,
                                clip=clip)

    colormap = {'log10': colors.LogNorm,
                'sqrt': sqrtwrapper,
                'none': colors.Normalize}

    # Are these correct?
    efilts = {'w': '',
              'b': '[energy=500:7000]',
              'u': '[energy=300:500]',
              's': '[energy=500:1200]',
              'm': '[energy=1200:2000]',
              'h': '[energy=2000:7000]'}

    # Work out the event file + VFS stacks to use before reading
    # anything else, since this lets us cache the result, which may
    # help out with really-large datasets (but only if the
    # band matches).
    #
    # TODO: look at this
    # Note that this cache is fairly pointless as we repeat this
    # code several times (for different scalings), when we should
    # perhaps just change the scaling here (to avoid re-creating
    # everything).
    #
    # key = (stack, cpt), value = file name incl VFS
    evt_name = {}

    # keys = file name incl VFS, values = number of times
    evt_count = defaultdict(int)

    for hull_idx, shull in enumerate(stackhulls):

        stack = shull['stack']
        cpt = shull['component']
        key = stack, cpt

        evtfile = evtfiles[stack]

        # What band to use?
        #
        bname = hulldata[key]['eband']

        if bname == 'w':
            bspec = '[bin sky=::32]'
        else:
            bspec = '[bin sky=::8]'

        efilt = efilts[bname]
        iname = evtfile + efilt + sfilt + bspec

        evt_name[key] = iname
        evt_count[iname] += 1

    evt_cache = {}

    for hull_idx, shull in enumerate(stackhulls):

        stack = shull['stack']
        cpt = shull['component']
        key = stack, cpt

        if hull_idx % nplots == 0:
            if page_idx > 0:
                save_plot(page_idx)

            page_idx += 1

            if page_idx == npages:
                nplots_in_page = nhulls - nplots * (page_idx - 1)
            else:
                nplots_in_page = nplots

            fig = plt.figure(figsize=(8, 8))

            # the only plot which we'll show axes for; it is intended
            # to be the bottom-left plot of the display.
            #
            nrows = int((nplots_in_page - 1) // nsize)
            axplot = nrows * nsize + 1

        plot_idx = 1 + (hull_idx % nplots)

        # pick the same scale when log scale is used
        # maybe...
        # if evtscale == 'log10':
        #     iopts['threshold'] = [0, 2]
        #

        # What band to use?
        #
        bname = hulldata[key]['eband']

        # Read in the date or grab it from the cache
        # *UNLESS* it is a special case.
        #
        if stack == 'acisfJ0534316p220052_001':
            # The > 4Gb Crab event file.

            ax = fig.add_subplot(nsize, nsize, plot_idx)
            ax.text(0.5, 0.8, stack, ha='center')
            ax.text(0.5, 0.3, 'No image as evt3 too large',
                    ha='center')
            continue

        iname = evt_name[key]
        if iname in evt_cache:
            cr = evt_cache[iname]

            evt_count[key] -= 1
            if evt_count[key] < 1:
                del evt_cache[iname]

        else:
            try:
                cr = pycrates.read_file(iname)
            except IOError as e:
                print("Unable to read {}".format(iname))
                raise e

            # cache the result if we re-use it
            if evt_count[key] > 1:
                evt_cache[iname] = cr

        # could use the transfrom stored in the hullmap, but use
        # the one read in from the image for now.
        #
        # tr = hullmap[stack][0]['transform']
        tr = cr.get_transform('EQPOS')
        tr2 = cr.get_transform('SKY')
        assert tr2 is not None
        wcs = tr2wcs(tr, tr2)

        # Get the limits. Could use all images in the region,
        # but this can be biased by counts in an unrelated
        # source (have seen this).
        #
        # So try and filter to just the counts within the hull,
        # using the SKY coordinate system.
        #
        hullbounds = hulldata[key]['pos']
        pixvals = cr.get_image().values
        ny, nx = pixvals.shape

        # create the logical axis; center of bottom-left pixel is
        # 1.0,1.0
        #
        # convert to sky coordinates
        #
        ly, lx = np.mgrid[1:ny + 1, 1:nx + 1]
        cs_log = 1.0 * np.vstack((lx.flatten(), ly.flatten()))
        cs_sky = tr2.apply([cs_log])[0]
        sx = cs_sky[0]
        sy = cs_sky[1]

        # Ideally would restrict to just the counts within the
        # hull (in case there's an unrelated bright source in the
        # field of view). In this case may want to ensure vmin is 0
        # just in csae the bounded hull has no low-value pixels in
        # it (seems unlikely, but a reasonable precaution).
        #
        fpixvals = filter_image(pixvals.flatten(), sx, sy, hullbounds)

        # One reason for not using the data to get the minimum is that
        # it is possible for pixel values to overflow, which results
        # in large negative values (unless we tell the DM to use a
        # larger int for the pixel values, but we already have memory
        # issues for some cases), and we do not want to include these
        # in the scaling.
        #
        # vmin = np.nanmin(fpixvals)
        vmin = 0
        vmax = np.nanmax(fpixvals)

        # It is possible for the filter to exclude all pixels if the
        # convex hull is small enough (given the binning of the pixels)
        # This has been seen with ens0000100_001.002, so in this case
        # ignore the filter as a work around.
        #
        if vmax == 0:
            vmax = np.nanmax(pixvals)

        if evtscale == 'log10' and vmin <= 0:
            # what to do with 0 elements?
            # Assumption is that the images are count arrays
            #
            vmin = 1

        cfunc = colormap[evtscale]
        cscale = cfunc(vmin=vmin, vmax=vmax)

        # Set up the plot data; note this uses the stack-specific
        # transform.
        #
        ax = fig.add_subplot(nsize, nsize, plot_idx, projection=wcs)
        ax.set_aspect('equal')
        ax_trans = ax.get_transform('world')

        if plot_idx == 1:
            ptitle = title + '  Page {}/{}'.format(page_idx,
                                                   npages)

            fig.text(0.5, 0.96, ptitle, fontsize=16,
                     color=tcol,
                     horizontalalignment='center')

            if qatitle is not None:
                fig.text(0.5, 0.92, qatitle, fontsize=16,
                         color=qacol,
                         horizontalalignment='center')

            if evtscale == 'none':
                sclbl = 'linear'
            else:
                sclbl = evtscale
            fig.text(0.95, 0.05, "scale={}".format(sclbl),
                     fontsize=lblsize,
                     horizontalalignment='right')

        # This happened during development (should since have been
        # fixed), but left in just in case.
        #
        if vmax <= vmin:
            print("INTERNAL ERROR vmin={} vmax={}".format(vmin, vmax) +
                  " file={}".format(iname))
            plt.text(0.5, 0.5, "SKIPPING IMAGE",
                     horizontalalignment="center",
                     transform=ax.transAxes)
            continue

        # for log scaling can have values <= 0; these turn into
        # transparent pixels, which is slightly distracting, so replace
        # them with a minimum value. As the images are count images,
        # then the minimum value should be log10(1) = 0
        #
        ivals = cr.get_image().values

        if np.all(ivals < 0):
            # have seen at least one strange case; add some hacky
            # debugging
            #
            print("WARNING: all values are negative... " +
                  "for file={}".format(iname))

        if evtscale == 'log10':
            # Should this really be ivals <= vmin as the condition?
            # The input image is assumed to be a counts image (i.e.
            # integer values), and we have vmin >= 1 for log10 scaling
            # from above.
            #
            ivals[ivals <= 0] = vmin

        # probably don't need vmin/vmax arguments, but as we have
        # them (they were added whilst debugging a problem that turned
        # out not to be a problem).
        #
        # Hide any warnings about square root of a negative number.
        # The colorbar should be safe from these but moved in during
        # testing (before realising that also needed to warn when
        # creating the hardcopy version) and have decided to leave
        # here.
        #
        oldsettings = np.seterr(all='ignore')
        try:
            plt.imshow(ivals, origin='lower', norm=cscale,
                       vmin=vmin, vmax=vmax)
            if colorbar:
                plt.colorbar()

        finally:
            np.seterr(**oldsettings)

        # Turn off autoscaling, so that we can draw on other hulls
        # and not worry about them overlapping anything.
        #
        ax.autoscale(enable=False)

        # Draw on any PSFs. Do this first so they appear under
        # everything else.
        #
        draw_xmdat3(xmdat3map[stack], ax, ax_trans)

        label_plot(stack, cpt, bname, ax,
                   hullmap, ensemblemap)

        if show_other_stack_hulls:
            add_hulls_from_other_stacks(stack, stacks, hullmap, ax)

        add_hulls(stack, cpt, hullmap,
                  master_hull, qahulls, ax,
                  master_color=master_color,
                  qa_color=qa_color)

        # clean up the plot.
        #
        labelize_plot(ax, plot_idx == axplot)

    # Don't forget to save the last page.
    save_plot(page_idx)
    assert page_idx == npages
    return {'npages': npages}


def draw_ensemble_outline(ensemble, mhulls, hulls, qas, fov3files,
                          qa_color='cyan'):
    """Draw up the FOVs for the ensembles and the different stack-level hulls.

    Parameters
    ----------
    ensemble : str
        The ensemble name; it is used for the plot title and error
        messages.
    mhulls : dict
        The master hulls - the second argument of
        utils.read_master_hulls - which is a dictionary keyed by
        the master id and the values give the master-hull details.
    hulls : list of list of stack hulls
        The stack-level hulls; each entry in the list is the output of
        read_hulls_from_mrgsrc3
    qas : None or dict
        If not None then the QA hulls for this ensemble, where the
        key is the master id and the value is the return value
        of chs_create_review_images_mpl:read_qa_hulls.
    fov3files : list of str
        The FOV3 files for the stacks in the ensemble. This is assumed
        to be the FOV files for those stacks containing hulls (i.e.
        they act as a reasonable bounding box for the display).
    qa_color : str, optional
        The color for any QA hulls and labels.

    """

    # Use the first transform we find as the base.
    #
    tr = None
    for hlist in hulls:
        for hull in hlist:
            tr = hull['transform']
            break

    if tr is None:
        raise IOError("No hulls found for " +
                      "ensemble: {}".format(ensemble))

    wcs = tr2wcs(tr)
    fig = plt.figure(figsize=(10, 10))
    ax = fig.add_subplot(1, 1, 1, projection=wcs)
    ax_trans = ax.get_transform('world')

    # We want a different color for each stack, but note that there
    # are multiple polygons per stack, so we can not just rely on
    # matplotlib's automatic cycling of colors.
    #
    # Use the default matplotlib set of colors for cycling.
    #
    colors = plt.rcParams['axes.prop_cycle'].by_key()['color']
    _nextcolor = cycler(color=colors)()

    def nextcolor():
        # return _nextcolor.next()['color']
        return next(_nextcolor)['color']

    for fovfile in fov3files:
        # assume everything is a polygon
        cr = pycrates.read_file(fovfile)
        color = nextcolor()
        for eqpos in cr.get_column('EQPOS').values:
            ra = eqpos[0]
            dec = eqpos[1]
            idx = np.isfinite(ra)

            ra = ra[idx]
            dec = dec[idx]

            plt.plot(ra, dec, color=color, alpha=0.07,
                     transform=ax_trans)

    # Draw these somewhat transparent so that the numbering added
    # later can be seen (it can get busy).
    #
    # stack_color = 'orange'
    # main_color = 'gold'
    stack_color = 'black'
    main_color = 'black'

    for i, hlist in enumerate(hulls):

        for stkhull in hlist:
            if stkhull['mancode']:
                lstyle = 'dotted'
            else:
                lstyle = 'solid'

            # should this use chs_status.DELETE?
            if stkhull['match_type'] == "deleted":
                use_color = 'red'
            else:
                use_color = stack_color

            eqpos = stkhull['eqpos']
            ra = eqpos[0]
            dec = eqpos[1]

            plt.plot(ra, dec, linewidth=1, linestyle=lstyle,
                     alpha=0.6,
                     transform=ax_trans,
                     color=use_color)

    # Now the master hulls; these are going to overwrite the stack-level
    # hulls for the majority of cases (single stack hull), so use a
    # thinner line width, which is not ideal.
    #
    title_col = 'k'
    for mid, mhull in mhulls.items():
        if chs_status.is_qa(mhull['status']):
            title_col = 'r'
            continue

        eqpos = mhull['eqpos']
        ra = eqpos[0]
        dec = eqpos[1]

        plt.plot(ra, dec, linewidth=1, linestyle='solid',
                 alpha=0.6,
                 transform=ax_trans,
                 color=main_color)

        label_hull(ax_trans, mhull,
                   "{}".format(mid),
                   main_color)

    # Do we have any qa hulls?
    #
    if qas is not None:
        for mid, qahulls in qas.items():

            # We label each "qa component"
            #
            nqas = len(qahulls)
            for i, qahull in enumerate(qahulls):
                draw_hull(ax_trans, qahull, qa_color,
                          1, 'dashed')

                lbl = "QA {}".format(mid)
                if nqas > 1:
                    lbl += " [{}]".format(i + 1)

                label_hull(ax_trans, qahull, lbl,
                           qa_color)

    # Report on the number of stacks with hulls, not the number
    # in the ensemble.
    #
    lbl = "Ensemble {}   #stacks={}  #hulls={}".format(ensemble,
                                                       len(fov3files),
                                                       len(mhulls))
    plt.title(lbl, fontsize=18, color=title_col)

    ax.coords['ra'].set_major_formatter('hh:mm:ss')
    ax.coords['dec'].set_major_formatter('dd:mm:ss')

    ax.set_aspect('equal')
