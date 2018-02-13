"""
Create review plots.

This requires CIAO
"""

import os

import numpy as np

import six

import pycrates

from pychips import all as pychips

from crates_contrib.utils import scale_image_crate

# from . import chs_utils as utils
import chs_utils as utils


def draw_hulls_and_images(master_hull,
                          stackhulls,
                          hullmap,
                          stkevt3dir,
                          outdir,
                          ensemble,
                          revision,
                          qahulls=None,
                          evtscale='log10',
                          master_color='gold',
                          qa_color='cyan',
                          axscale=0.5,
                          showhack=False):
    """Draw individidual hull+image for all members of an ensemble.

    Parameters
    ----------
    master_hull: dict
        Contains the master hull: fields are 'master_id', 'status',
        and 'eqpos'. The 'status' field should be either 'qa' or
        'okay'.
    stackhulls : list
        What stack-level hulls form this master hull? Each entry
        is a dictionary with the keys 'stack' and 'component'.
    hullmap : dict
        The stack-level hull data, stored by the stack id. The values
        are those returned by read_hulls_from_mrgsrc3.
    stkevt3dir : str
        The location of the stkevt3 files.
    outdir : str
        The output directory, which must exist.
    ensemble : str
        The ensemble value.
    revision : int
        The revision number of the file.
    qahulls : None or list of dict, optional
        This is only used if master_hull['status'] is set to 'qa'.
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
    showhack : bool, optional
        Does the window have to be displayed before it can be printed?
        This seems to depend on the system.

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

    hulldata = {}
    for stack, hulls in six.iteritems(hullmap):
        if stack not in stacks:
            continue

        for hull in hulls:
            assert hull['stack'] == stack
            key = stack, hull['component']
            assert key not in hulldata
            hulldata[key] = hull

    nhulls = len(hulldata)

    # Need limits for the plots: it would be easier if we could just
    # rely on ChIPS to calculate this, but there's no guarntee that
    # the limits will include all the data once the aspect ratio is
    # set.
    #
    # Technically the QA and master hulls should be included here, but
    # they should all, by definition, be no larger than the
    # stack-level hulls.
    #
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

    # pick the "middle" of the shape. need to calculate the
    # zoom factor per plot, which is less than ideal.
    #
    ra_lims = np.asarray(ra_lims)
    dec_lims = np.asarray(dec_lims)

    # dra = ra_lims.ptp()
    # ddec = dec_lims.ptp()

    # ra_center = ra_lims.min() + dra / 2
    # dec_center = dec_lims.min() + ddec / 2

    ra_min = ra_lims.min()
    ra_max = ra_lims.max()
    dec_min = dec_lims.min()
    dec_max = dec_lims.max()

    dra = ra_max - ra_min
    ddec = dec_max - dec_min

    # expand limits on each side
    #
    ra_min -= axscale * dra
    ra_max += axscale * dra

    dec_min -= axscale * ddec
    dec_max += axscale * ddec

    dra = ra_max - ra_min
    ddec = dec_max - dec_min
    dec_mid = (dec_min + dec_max) / 2.0
    dra = dra * np.cos(dec_mid * np.pi / 180.0)

    # With the current set up, the image data is fixed to the limits
    # that are being shown, so zooming out doesn't really help.
    #
    # zoomfactor = 0.9
    zoomfactor = 1.0

    if ddec > dra:
        def reset_limits():
            pychips.limits(pychips.Y_AXIS, dec_min, dec_max)
            pychips.zoom(zoomfactor)

    else:
        def reset_limits():
            pychips.limits(pychips.X_AXIS, ra_min, ra_max)
            pychips.zoom(zoomfactor)

    # It's not worth reading in the whole event file/image if
    # possible.
    #
    sfilt = utils.make_spatial_filter_range(ra_min, ra_max,
                                            dec_min, dec_max)

    # Find the event files for the stacks.
    #
    evtfiles = {}
    for stack in stacks:
        path = os.path.join(stkevt3dir,
                            '{}*_evt3.fits*'.format(stack))
        evtfiles[stack] = utils.find_single_match(path)

    # Iterate over the stack, so that multiple hulls from the
    # same stack are shown in the same plot.
    #
    # TODO: is this true?
    #
    nstacks = len(stacks)

    # could go for a more-adaptive scheme
    nsize = np.ceil(np.sqrt(nstacks)).astype(np.int)
    nplots = nsize * nsize
    if nsize > 3:
        nsize = 3
        nplots = nsize * nsize
        npages = 1 + nstacks // nplots
    else:
        npages = 1

    def save_plot(pagenum):
        """Create PNG output then destroy the window"""

        # It looks like we need to set the display to True to get
        # any sensible hardcopy output.
        #
        # pychips.set_window(['display', True])  is the comment still true?
        fmt = 'hull.{}.{:03d}.p{:03d}.v{:03d}.{}.png'
        outfile = fmt.format(ensemble, masterid, pagenum, revision,
                             evtscale)

        outfile = os.path.join(outdir, outfile)
        if showhack:
            pychips.set_window(['display', True])

        pychips.print_window(outfile, ['clobber', True])
        print("Created: {}".format(outfile))
        pychips.delete_window()

    title = "{} {:03d}".format(ensemble, masterid) + \
        "  #stacks={} #hulls= {}".format(nstacks, nhulls)

    if master_hull['status'] != 'okay':
        title = "{{\color{{red}}{}}} {}".format(master_hull['status'].upper(),
                                                title)

    # should not be needed, but appears to be
    hulldepth = 200

    # TODO: should we draw on the FOV?
    def draw_hull(hull, hcolor, lthick, lstyle):
        "line-style handling is messy"

        if 'mancode' in hull:
            if hull['mancode'] == 0:
                lstyle = 'solid'
            else:
                lstyle = 'dot'

        ropts = {'*.color': hcolor,
                 'opacity': 1,
                 'edge.style': lstyle,
                 'edge.thickness': lthick,
                 'depth': hulldepth,
                 'fill.style': 'nofill'}
        copts = {'*.color': hcolor,
                 'line.style': lstyle,
                 'line.thickness': lthick,
                 'depth': hulldepth,
                 'symbol.style': 'none'}

        # the regions may be degenerate, so fall through to a
        # curve if this happens. Should really only catch that
        # case, to avoid hiding other problems, but leave that
        # for now.
        #
        eqsrc = hull['eqpos']
        x = eqsrc[0]
        y = eqsrc[1]
        try:
            pychips.add_region(x, y, ropts)
        except RuntimeError:
            pychips.add_curve(x, y, copts)

    page_idx = 0
    nplots_in_page = None
    for stack_idx, stack in enumerate(stacks):

        if stack_idx % nplots == 0:
            if page_idx > 0:
                save_plot(page_idx)

            page_idx += 1

            if page_idx == npages:
                nplots_in_page = nstacks - nplots * (page_idx - 1)
            else:
                nplots_in_page = nplots

            pychips.add_window(8, 8, 'inches', ['display', False])
            pychips.add_frame()

            for i in range(nplots_in_page):
                # useful to have hard-coded names; even if this is the
                # default value
                pychips.add_plot(['id', 'plot{}'.format(i + 1)])

                """ Can not create the plots with the data we have

                WHY IS THIS?

                # TODO: do the ra limits need reversing?
                pychips.add_axis(pychips.XY_AXIS, 0,
                                 ra_min, ra_max, dec_min, dec_max, tr)
                pychips.set_data_aspect_ratio('1:1')
                pychips.set_axis(['automin', False, 'automax', False,
                                  'tickstyle', 'outside'])
                pychips.zoom(zoom)
                """

            # the only plot which we'll show axes for; it is intended to
            # be the bottom-left plot of the display.
            #
            nrows = int((nplots_in_page - 1) // nsize)
            axplot = nrows * nsize + 1

            pychips.grid_objects(nsize, nsize, 0, 0, 0)

            # add a title to the first plot
            pychips.set_current_plot('plot1')
            ptitle = title + '  Page {}/{}'.format(page_idx,
                                                   npages)

            ptitle = ptitle.replace('_', '\\_')
            pychips.add_label(0.5, 0.92, ptitle,
                              {'coordsys': pychips.FRAME_NORM,
                               'halign': 0.5,
                               'size': 24})

        plot_idx = 1 + (stack_idx % nplots)
        pychips.set_current_plot('plot{}'.format(plot_idx))

        # Set up the plot data; note this uses the stack-specific
        # transform.
        tr = hullmap[stack][0]['transform']
        pychips.add_axis(pychips.XY_AXIS, 0,
                         ra_min, ra_max, dec_min, dec_max, tr)
        pychips.set_data_aspect_ratio('1:1')
        pychips.set_axis(['automin', False, 'automax', False,
                          'tickstyle', 'outside'])

        # draw on the image, below everything
        #
        iopts = {'depth': 50}

        # pick the same scale when log scale is used
        # maybe...
        # if evtscale == 'log10':
        #     iopts['threshold'] = [0, 2]
        #
        evtfile = evtfiles[stack]

        # What band to use?
        #
        # TODO: what to do when multiple hulls per stack?
        #       in fact, do we get all hulls from this stack in a
        #       single plot?
        #
        bnames = set([])
        for hulls in hullmap[stack]:
            bnames.add(hulls['band'])

        if len(bnames) == 1:
            bname = bnames.pop()
        else:
            raise RuntimeError("BANDS={}".format(bnames))

        if bname == 'w':
            bspec = '[bin sky=::64]'
        else:
            bspec = '[bin sky=::8]'

        # Are these correct?
        efilts = {'w': '',
                  'b': '[energy=500:7000]',
                  'u': '[energy=300:500]',
                  's': '[energy=500:1200]',
                  'm': '[energy=1200:2000]',
                  'h': '[energy=2000:7000]'}

        efilt = efilts[bname]
        iname = evtfile + efilt + sfilt + bspec
        try:
            cr = pycrates.read_file(iname)
        except IOError as e:
            print("Unable to read {}".format(iname))
            raise e

        oldsettings = np.seterr(all='ignore')
        try:
            scale_image_crate(cr, evtscale)
        finally:
            np.seterr(**oldsettings)

        pychips.add_image(cr, iopts)

        pychips.add_label(0.05, 0.9, "band={}".format(bname),
                          {'color': 'white', 'size': 20,
                           'coordsys': pychips.PLOT_NORM})

        pychips.add_label(0.95, 0.9, evtscale,
                          {'color': 'white', 'size': 20,
                           'halign': 1.0,
                           'coordsys': pychips.PLOT_NORM})

        # THIS SHOULD NOT BE NEEDED, BUT SEEMS TO BE FOR ONE CASE
        #
        pychips.shuffle_image(pychips.chips_back)

        # draw the hulls for this stack
        #
        # draw the other hulls first, as reference
        # (i.e. those from other stacks)
        #
        # NOTE: this draws on all hulls, so can be useful
        #       if nearby ones overlap.
        #
        for ostack in stacks:
            if ostack == stack:
                continue

            for hull in hullmap[ostack]:
                draw_hull(hull, 'steelblue', 1, 'dot')

        for hull in hullmap[stack]:
            draw_hull(hull, 'orange', 3, 'solid')

        # do we have a master hull to add?
        #
        if master_hull['status'] == 'okay':
            draw_hull(master_hull, master_color, 3, 'solid')
        else:
            for qahull in qahulls:
                draw_hull(qahull, qa_color, 3, 'longdash')

        # clean up the plot.
        #
        if plot_idx == axplot:
            pychips.set_xaxis(['tickformat', '%ra'])
            pychips.set_yaxis(['tickformat', '%dec'])
        else:
            pychips.hide_axis()

        reset_limits()

    # need to save the last page
    save_plot(page_idx)


_color = ['gray', 'cyan', 'steelblue', 'brown', 'beige']


def nextcolor(coldict):
    """Return a "new" color (only a small number are available).

    Parameters
    ----------
    coldict : dict
        The first time called this should be an empty dictionary.
        Pass in the same dictionary to the next call.

    """

    try:
        colnum = coldict['counter']
    except:
        colnum = 0

    try:
        color = _color[colnum]
    except IndexError:
        colnum = 0
        color = _color[colnum]

    coldict['counter'] = colnum + 1
    return color


def draw_ensemble_outline(ensemble, hulls, fov3files):
    """Draw up the FOVs for the ensembles and the different hulls.

    This creates a new ChIPS window but does not set the display
    to True.

    Parameters
    ----------
    ensemble : str
        The ensemble name; it is used for the plot title and error
        messages.
    hulls : list of list of stack hulls
        The stack-level hulls; each entry in the list is the output of
        read_hulls_from_mrgsrc3
    fov3files : list of str
        The FOV3 files for the stacks in the ensemble. This is assumed
        to be the FOV files for those stacks containing hulls (i.e.
        they act as a reasonable bounding box for the display).

    """

    pychips.add_window(10, 10, 'inches', {'display': False})

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

    pychips.add_axis(pychips.XY_AXIS, 0, 0, 1, tr)

    # Draw the hulls first (using a curve) so get auto-scaling,
    # with a high depth so they are above the FOVs, which will get
    # added later.
    #
    copts = {'depth': 150, 'symbol.style': 'none',
             'line.thickness': 2, '*.color': 'orange'}
    nhulls = 0
    for hlist in hulls:
        for stkhull in hlist:
            if stkhull['mancode']:
                copts['line.style'] = 'dot'
            else:
                copts['line.style'] = 'solid'

            eqpos = stkhull['eqpos']
            ra = eqpos[0]
            dec = eqpos[1]
            pychips.add_curve(ra, dec, copts)
            nhulls += 1

    """
    # Note that the axis limits are chosen based on the hull location,
    # not the extent of the FOV files. There are several ways to do
    # this (e.g. draw the hulls as regions), but for now I'm moving
    # the axis limits *before* drawing the FOV files.
    #
    xr1 = pychips.get_plot_xrange()
    yr1 = pychips.get_plot_yrange()
    pychips.set_data_aspect_ratio('1:1')

    ymid = (yr1[0] + yr1[1]) / 2.0
    xscale = np.cos(ymid * np.pi / 180.0)

    dx1 = (xr1[0] - xr1[1]) * xscale
    dy1 = yr1[1] - yr1[0]
    if dy1 > dx1:
        pychips.limits(pychips.Y_AXIS, yr1[0], yr1[1])
    else:
        pychips.limits(pychips.X_AXIS, xr1[0], xr1[1])
    """

    # Use a region for the FOV file so that the opacity of the
    # field edges can be reduced.
    #
    ropts = {'depth': 50, 'fill.style': 'none', 'opacity': 0.2}

    colinfo = {}
    ropts['*.color'] = nextcolor(colinfo)

    xr1 = pychips.get_plot_xrange()
    yr1 = pychips.get_plot_yrange()

    for fovfile in fov3files:
        # assume everything is a polygon
        cr = pycrates.read_file(fovfile)
        for eqpos in cr.get_column('EQPOS').values:
            ra = eqpos[0]
            dec = eqpos[1]
            idx = np.isfinite(ra)

            ra = ra[idx]
            dec = dec[idx]

            # Hopefully there's no case where this call fails
            # with FOV files.
            pychips.add_region(ra, dec, ropts)

            # xr1 is swapped
            rmax = ra.max()
            rmin = ra.min()
            if rmin < xr1[1]:
                xr1[1] = rmin
            if rmax > xr1[0]:
                xr1[0] = rmax

            dmax = dec.max()
            dmin = dec.min()
            if dmin < yr1[0]:
                yr1[0] = dmin
            if dmax > yr1[1]:
                yr1[1] = dmax

        ropts['*.color'] = nextcolor(colinfo)

    # Try calculating the boundaries after adding in the regions,
    # although now have to rely on manually calculating the desired
    # range since regions do not cause the axes to change size.
    #
    pychips.set_data_aspect_ratio('1:1')

    ymid = (yr1[0] + yr1[1]) / 2.0
    xscale = np.cos(ymid * np.pi / 180.0)

    dx1 = (xr1[0] - xr1[1]) * xscale
    dy1 = yr1[1] - yr1[0]
    if dy1 > dx1:
        pychips.limits(pychips.Y_AXIS, yr1[0], yr1[1])
    else:
        pychips.limits(pychips.X_AXIS, xr1[0], xr1[1])

    # Use the same title size as used in draw_hulls_and_images
    nstacks = len(hulls)
    lbl = "Ensemble {}   #stacks={}  #hulls={}".format(ensemble,
                                                       nstacks,
                                                       nhulls)
    lbl = lbl.replace('_', '\\_')
    pychips.set_plot_title(lbl)
    pychips.set_plot({'title.size': 24})
    pychips.set_xaxis({'tickformat': '%ra'})
    pychips.set_yaxis({'tickformat': '%dec'})

    # pychips.set_window({'display': True})
