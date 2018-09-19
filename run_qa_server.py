#!/usr/bin/env python

"""

Usage:

  ./run_qa_server.py datadir [port]

The default port is 8070. The datadir is the location of the
master-hull data (the location of the ens<> directories).

The webserver assets are found relative to the location of the
script.

Data will be written to the working directory.

Things to do

 - hull review page: regions could be added to JS when the page is
   created, which would avoid some of the callbacks needed
 - should draw all stack hulls in JS9 window

  These two are almost done; need to remove or change the behavior of the
  reload region button, then can remove/clean up the region endpoint
  and need to send in the hull centers for the master hulls

  Need to start labelling the master hulls (for when there are qa cases)
  so that can delete a single hull.

  How to delete the convex hull after adjustment? Maybe don't do this
  ambiently, but wait for the user to save the hull, in which case the
  display can be updated at that time (or perform some other action).

 - when the master hull is changed, should change the JS data
   representation and then use that to update other JS9 displays.

 - why is ens0008300_001.013 pre-labelled as manual QA?
   (may only appear on re-loads on firefox?)
 - do we want to show all masters and other hulls from the stack
   in the JS9 window? Now add in all the hulls from the stack.

 - are HRC events being binned correctly (in particular, image size)?

 - handle band selection (hard for case where multiple components)
    - want easy band choice from display too
    - not sure what the best UI is, particularly if end up with
      "combined" bands

 - implement any API needed for SAMP (low priority for now)
 - implement read/save text and decisions
 - send in suggested low/high pixel value? Or calculate it directly

 - UI to hide/show the different layers *(e.g. if include PSF regions)

Notes:

This was based on the skeleton code provided at

  https://gist.github.com/bradmontgomery/2219997

See also https://pymotw.com/2/BaseHTTPServer/

"""

import glob
import gzip
import json
import os
import sys

from six.moves.BaseHTTPServer import BaseHTTPRequestHandler, \
    HTTPServer

from six.moves.urllib.parse import urlparse

import numpy as np

from jinja2 import Environment, FileSystemLoader, select_autoescape, \
    TemplateNotFound

try:
    import pycrates
except ImportError:
    sys.stderr.write("Needs crates!\n")
    sys.exit(1)

import chs_utils as utils


class CHSServer(HTTPServer):
    """Allow information to be passed to the handler.

    The store argument should be a dictionary with
    keywords:
       datadir  - the directory containing the ensemble/master hull
                  information (this is created by the CHS
                  master-match pipeline)
       evt3dir  - the location of the stack event files
       userdir  - the location for storing user information
       xmdatdir - the location of the xmdat files (stored in subdirectories
                  labelled by the stack name).
       webdir   - the location of the web assets
       environment - the Jinja2 environment

    These must all be existing directories.
    """

    def __init__(self, store, *args, **kwargs):
        context = {}
        for k in ['datadir', 'evt3dir', 'userdir', 'xmdatdir',
                  'webdir']:
            try:
                d = store[k]
            except KeyError:
                raise ValueError("Missing key={} in store".format(k))

            if not os.path.isdir(d):
                raise ValueError("Not a directory: {}".format(d))

            context[k] = d

        context['environment'] = store['environment']

        HTTPServer.__init__(self, *args, **kwargs)
        self.context = context


def valid_ensemble(datadir, ensemble):
    """Should we treat this as an ensemble?

    A valid ensemble value is defined as one in which the name
    matches ens???????_??? and the directory exists (where the
    question marks represent 0-9.

    Parameters
    ----------
    datadir : str
        The path to the directory containing the ensemble-level
        products.
    ensemble : str
        The ensemble id.

    Returns
    -------
    valid : bool

    """

    if len(ensemble) != 14 or not ensemble.startswith('ens') or \
            ensemble[-4] != '_':
        return False

    # Ensure that in ensXXXXXXX_YYY parts, XXX and YYY are integers
    try:
        int(ensemble[3:-4])
        int(ensemble[-3:])
    except ValueError:
        return False

    path = os.path.join(datadir, ensemble)
    return os.path.isdir(path)


def parse_datadir(datadir, userdir):
    """Extract useful information from the diretory.

    Parameters
    ----------
    datadir : str
        The path to the directory containing the ensemble-level
        products.
    userdir : str
        The path to the directory containing the user's decisions.

    Returns
    -------
    state : dict or None
        A dictionary or None, if no data can be found there.

    Notes
    -----
    At present the contents of datadir are highly structured;
    let's see if this will need relaxing.
    """

    pat = os.path.join(datadir, "ens*")
    ensembledirs = glob.glob(pat)
    if len(ensembledirs) == 0:
        utils.errlog("no ensemble dirs found in {}".format(datadir))
        return None

    """
    # HACK: use a subset of ensembles to save time
    #
    nmax = 60
    dbg("restricting to first {} ensembles".format(nmax))
    ensembledirs = ensembledirs[:nmax]
    """

    # Get the latest version of each ensemble
    #
    store = {}
    for ensembledir in ensembledirs:
        ensemble = ensembledir.split('/')[-1]
        # UGH: this is ugly; note the /../ added to the path
        # TODO: fix this
        jcts = utils.read_ensemble_json(ensembledir + "/../",
                                        userdir,
                                        ensemble)
        if jcts is None:
            store[ensemble] = {'error': True, 'ensemble': ensemble}
            continue

        v = jcts['latest_version']
        info = jcts['versions'][v]
        e = info['name']
        store[e] = info

    # return in ensemble order
    keys = sorted(store.keys())
    return [store[k] for k in keys]


def get_data_summary(datadir, userdir):
    """Return the data about all the available ensembles.

    Parameters
    ----------
    datadir : str
        The path to the directory containing the ensemble-level
        products.
    userdir : str
        The path to the directory containing the user's decisions.

    Returns
    -------
    summary : dict or None
    """

    state = parse_datadir(datadir, userdir)
    if state is None:
        return None

    # split up the ensembles
    #
    # QUESTION: what are the valid status values?
    #
    out = {'todos': [],
           'reviews': [],
           'completed': [],
           # 'manual': [],
           'errors': [],
           'usernotes': '',
           'lastmodified': ''}

    for ens in state:
        if 'error' in ens:
            out['errors'].append(ens['ensemble'])
            continue

        status = utils.get_user_setting(ens, 'status')
        if status == 'completed':
            key = 'completed'
        elif status == 'review':
            key = 'reviews'
        else:
            key = 'todos'

        out[key].append(ens)

    # Has the user saved any notes?
    infile = os.path.join(userdir, 'summary.json')
    if os.path.isfile(infile):
        jcts = utils.read_json(infile)
        if jcts is not None:
            out['usernotes'] = jcts['usernotes']
            out['lastmodified'] = jcts['lastmodified']

    # How about the datatable settings?
    infile = os.path.join(userdir, 'datatable.json')
    if os.path.isfile(infile):
        jcts = utils.read_json(infile)
        if jcts is not None:
            out['datatable'] = jcts

    return out


def read_ds9_region(infile):
    """Read in DS9 region file and turn into JS9 string.

    This is specialized to the CHS polygon files, where only the first
    polygon shape is used. *LIMITED* validation.
    """

    try:
        with open(infile, 'r') as fh:
            csys = None
            for l in fh.readlines():
                l = l.strip()
                if l.startswith('#') or l.startswith('global '):
                    continue

                if csys is None:
                    if l.find(' ') != -1:
                        utils.errlog("expected coordsys, found " +
                                     "[{}] in {}".format(l, infile))
                        return None

                    csys = l
                    continue

                if not l.startswith('polygon('):
                    utils.errlog("expected polygon(...), found " +
                                 "[{}] in {}".format(l, infile))
                    return None

                return "{}; {}".format(csys, l)

    except IOError as exc:
        utils.errlog("can not read region file:" +
                     "\n{}\n{}".format(infile, exc))
        return None


def regstr_to_coords(regstr):
    """Extract RA,Dec vertices from a polygon string.

    Parameters
    ----------
    regstr : str
        The output of read_ds9_region. Can not be None.

    Returns
    -------
    ra, dec, label : list of real, list of real, str or None
        The RA and Dec values from the polygon in regstr,
        and any label in the regstr. For now label is always None

    Raises
    ------
    ValueError
        When regstr does not match the very-limited supported format
        of 'fk5; polygon(...)...# text={...} ...'.

    Notes
    -----
    Really we should be reading the coordinates dircectly from a file,
    but at present we don't have the file to hand.
    """

    if regstr is None:
        raise ValueError("regstr is None")
    if not regstr.startswith('fk5; polygon('):
        raise ValueError("expected 'fk5; polygon(...' but " +
                         "got '{}'".format(regstr))

    epos = regstr.find(')')
    if epos == -1:
        raise ValueError("No end ) in '{}'".format(regstr))

    # Coordinates are labelled with a trailing d.
    def clean(c):
        if c.endswith('d'):
            v = c[:-1]
        else:
            v = c

        try:
            return float(v)
        except ValueError:
            raise ValueError("Expected number+d but got '{}'".format(c))

    cs = regstr[13:epos].split(',')
    npts = len(cs)
    if npts % 2 == 1:
        raise ValueError("Missing or extra coordinate value in '{}'".format(regstr))

    nv = npts // 2
    cs = np.asarray([clean(c) for c in cs]).reshape(nv, 2).T

    return list(cs[0]), list(cs[1]), None


def apply_template(env, tmplname, args):

    try:
        tmpl = env.get_template(tmplname)
    except TemplateNotFound as exc:
        utils.errlog("No template for " + tmplname)
        out = "<!DOCTYPE html><html><head><title>ERROR</title>"
        out += "</head><body><p>Unable to find template:"
        out += "<strong>{}</strong>\n".format(tmplname)
        out += "Error: {}".format(exc)
        out += "</body></html>"
        return 404, out

    out = tmpl.render(**args)
    return 200, out


def create_index_page(env, datadir):
    """Create the top-level page: ensemble status.

    Parameters
    ----------
    env
        Jinja2 environment
    datadir : str
        The path to the directory containing the ensemble-level
        products.

    Returns
    -------
    status, page : int, str
        The response status and HTML contents
    """

    return apply_template(env, 'index.html', {'datadir': datadir})


def create_ensemble_index_page(env, datadir, userdir, ensemble):
    """Create the top-level ensemble page.

    Parameters
    ----------
    datadir : str
        The path to the directory containing the ensemble-level
        products.
    userdir : str
        The path to the directory containing the user's decisions.
    ensemble : str
        The ensemble id.

    Returns
    -------
    status, page : int, str
        The response status and HTML contents
    """

    """
    I guess we could error out here rather than when the AJAX
    call fails, but let's see how this works

    state = read_ensemble_json(datadir, userdir, ensemble)
    if state is None:
        out = "<!DOCTYPE html><html><head><title>ERROR</title>"
        out += "</head><body><p>There was an error when processing "
        out += "the directory {} for ".format(datadir)
        out += "ensemble {}".format(ensemble)
        out += "</body></html>"
        return 404, out

    """

    # Could cache the username
    return apply_template(env, 'ensemble.html',
                          {'ensemble': ensemble,
                           'username': os.getlogin()})


def eqpos_to_dict(eqpos, nvertex):
    """Create the dictionary representing this polygon.

    Parameters
    ----------
    eqpos : NumPy array
        This is a 2 by n array giving the ra, dec values of
        the polygon.
    nvertex : int
        The number of vertexes for the polygon (assumed to
        close the polygon). It must be less than or equal
        to the size of the second axis of eqpos (i.e. n).

    Returns
    -------
    dict : dictionary
        The keys are 'ra' and 'dec' which contain Python
        lists (*not* NumPy arrays, as they can not be serialized
        to JSON), limited to nvertex points each, and
        'ra0', 'dec0' which given an estimate of the center
        of the hull (not intended to be very accurate).

    """

    assert nvertex <= eqpos.shape[1], \
        "{} and {}".format(nvertex, eqpos.shape)

    ra = eqpos[0, :nvertex]
    dec = eqpos[1, :nvertex]

    # for now go with a *very* simple estimate
    ra0 = 0.5 * (ra.min() + ra.max())
    dec0 = 0.5 * (dec.min() + dec.max())

    return {'ra': list(ra), 'dec': list(dec),
            'ra0': ra0, 'dec0': dec0}


def create_master_hull_page(env,
                            datadir, userdir, rawdir, xmdatdir,
                            ensemble, revision, masterid):
    """Create the review page for a master hull.

    This is not ideal, since the master id depends on the version;
    i.e. v001 may have 2 masters but v002 may have 1 or 3. In fact,
    the master id number is not very stable - it can end up referring
    to different master hulls in different versions.

    Parameters
    ----------
    env
        The Jinja2 environment.
    datadir : str
        The path to the directory containing the ensemble-level
        products.
    userdir : str
        The path to the directory containing the user's decisions.
    rawdir : str
        The path to the actual master hull data (fits and region files)
    xmdatdir : str
        The path to the xmdat directory (the files are stored in
        sub-directories labelled by the stack name).
    ensemble : str
        The ensemble id.
    revision : int
        The revision.
    masterid : int
        The master hull

    Returns
    -------
    status, page : int, str
        The response status and HTML contents

    Notes
    -----
    Given that the version is fixed, this is a mixture of filling in
    the HTML in Python and in Javascript.
    """

    # TODO: use chs_utils.read_master_hulls and convert this code

    mid = "{:03d}".format(masterid)
    revstr = "{:03d}".format(revision)

    state = utils.read_ensemble_json(datadir, userdir, ensemble)
    if state is None:
        out = "<!DOCTYPE html><html><head><title>ERROR</title>"
        out += "</head><body><p>There was an error when processing "
        out += "the directory {} for ".format(datadir)
        out += "ensemble {} hull {}".format(ensemble, masterid)
        out += "</body></html>"
        return 404, out

    # Is this the latest version?
    #
    is_latest = revstr == state['latest_version']

    # Get the hull.
    #
    info = state['versions'][revstr]
    hulls = [h for h in info['masters']
             if h['masterid'] == mid]
    if len(hulls) == 0:
        # may need a better message for when switching between versions
        out = "<!DOCTYPE html><html><head><title>ERROR</title>"
        out += "</head><body><p>Invalid master hull.</p></body></html>"
        return 404, out

    elif len(hulls) > 1:
        utils.log("hulls = {}".format(hulls))
        out = "<!DOCTYPE html><html><head><title>INTERNAL ERROR</title>"
        out += "</head><body><p>FOUND MULTIPLE COPIES - see Doug!</p></body></html>"
        return 404, out

    hull = hulls[0]

    # Read in the hull data for both all the masters and all the
    # stacks in this version.
    #
    # There is repeated work being done here which could/should be
    # rationalised.
    #
    dname = os.path.join(rawdir, ensemble)
    if not os.path.isdir(dname):
        utils.errlog("missing rawdir {}".format(dname))
        out = "<!DOCTYPE html><html><head><title>INTERNAL ERROR</title>"
        out += "</head><body><p>Missing dir={}".format(dname)
        out += "- see Doug!</p></body></html>"
        return 404, out

    hullfile = os.path.join(dname,
                            utils.make_mhull_name(ensemble, revision))
    if not os.path.isfile(hullfile):
        utils.errlog("missing hullfile {}".format(hullfile))
        out = "<!DOCTYPE html><html><head><title>INTERNAL ERROR</title>"
        out += "</head><body><p>Missing hullfile={}".format(hullfile)
        out += "- see Doug!</p></body></html>"
        return 404, out

    # From the hull file,
    #
    # *) HULLMATCH block
    #
    # find those stacks associated with a master source, and the
    # energy band.
    #
    # *) HULLLIST block
    #
    # find the state of each master hull and, where possible,
    # read in the vertices.
    #
    # The component values are corrected for the COMPZERO keyword
    # (if set).
    #
    try:
        ds = pycrates.CrateDataset(hullfile, mode='r')
    except IOError as exc:
        utils.errlog("unable to read hullfile {} - {}".format(hullfile, exc))
        out = "<!DOCTYPE html><html><head><title>INTERNAL ERROR</title>"
        out += "</head><body><p>Error reading "
        out += "hullfile={}\nreason=\n{}".format(hullfile, exc)
        out += "- see Doug!</p></body></html>"
        return 404, out

    try:
        cr = ds.get_crate('HULLMATCH')
    except IndexError as exc:
        utils.errlog("unable to read HULLMATCH from " +
                     "hullfile {} - {}".format(hullfile, exc))
        out = "<!DOCTYPE html><html><head><title>INTERNAL ERROR</title>"
        out += "</head><body><p>Error reading HULLMATCH block of "
        out += "hullfile={}\nreason=\n{}".format(hullfile, exc)
        out += "- see Doug!</p></body></html>"
        return 404, out

    # We have aleady got this information, but recreate it
    #
    nstks = cr.get_key_value('STKIDNUM')
    if nstks is None:
        utils.errlog("missing STKIDNUM keyword in hullfile " +
                     "{}".format(hullfile))
        out = "<!DOCTYPE html><html><head><title>INTERNAL ERROR</title>"
        out += "</head><body><p>No STKIDNUM keyword in "
        out += "hullfile={}\n".format(hullfile)
        out += "- see Doug!</p></body></html>"
        return 404, out

    stack_map = {}
    for i in range(nstks):
        key = "STKID{:03d}".format(i)
        val = cr.get_key_value(key)
        if val is None:
            utils.errlog("missing {} keyword in ".format(key) +
                         "hullfile {}".format(hullfile))
            out = "<!DOCTYPE html><html><head><title>INTERNAL ERROR</title>"
            out += "</head><body><p>No {} keyword in ".format(key)
            out += "hullfile={}\n".format(hullfile)
            out += "- see Doug!</p></body></html>"
            return 404, out

        stack_map[val] = i

    cpt0 = cr.get_key_value('COMPZERO')
    if cpt0 is None:
        cpt0 = 0

    # Do not need the eband info to be indexed by master id
    #
    detectors = set([])
    hull_store = {}
    mid_by_component = {}
    ebands_by_component = {}
    mancode_by_component = {}
    likelihood_by_component = {}
    svdqa_by_component = {}
    centroid_by_component = {}
    mrgrev_by_component = {}
    for vals in zip(cr.Master_Id.values,
                    cr.STACKID.values,
                    cr.COMPONENT.values,
                    cr.LIKELIHOOD.values,
                    cr.STKSVDQA.values,
                    cr.EBAND.values,
                    cr.INCLUDE_IN_CENTROID.values,
                    cr.MRG3REV.values,
                    cr.MAN_CODE.values):

        midval, stackid, cpt, lkhood, svdqa, eband, \
            centroid, mrgrev, mancode = vals
        cpt -= cpt0

        if stackid.startswith('acisf'):
            detectors.add('acis')
        elif stackid.startswith('hrcf'):
            detectors.add('hrc')
        else:
            utils.errlog("unexpected stackid={}".format(stackid))
            out = "<!DOCTYPE html><html><head><title>INTERNAL ERROR</title>"
            out += "</head><body><p>Unsupported stackid="
            out += "{}\n".format(stackid)
            out += "- see Doug!</p></body></html>"
            return 404, out

        # NOTE: mancode is encoded "strangely", so remove the
        #       extra dimension
        if mancode.shape == (1,):
            mancode = mancode[0]

        # mancode is also stored as np.uint8 which seems to be problematic
        # to serialize to JSON, so explicitly convert to a "normal" int
        #
        # there's something similar going on with mrg3rev, but this time
        # it appears to be being written out as np.int64.
        #
        mancode = int(mancode)
        mrgrev = int(mrgrev)

        # Note: dropping match type, which is probably okay
        try:
            store = hull_store[midval]
        except KeyError:
            store = {'components': {}}
            hull_store[midval] = store

        try:
            store['components'][stackid].append(cpt)
        except KeyError:
            store['components'][stackid] = [cpt]

        key = (stackid, cpt)
        for check in [mid_by_component, ebands_by_component,
                      mancode_by_component,
                      likelihood_by_component, svdqa_by_component,
                      centroid_by_component, mrgrev_by_component]:
            assert key not in check, "key={}".format(key)

        mid_by_component[key] = midval
        ebands_by_component[key] = eband
        mancode_by_component[key] = mancode
        likelihood_by_component[key] = lkhood
        svdqa_by_component[key] = svdqa
        centroid_by_component[key] = centroid
        mrgrev_by_component[key] = mrgrev

    cr = None

    try:
        cr = ds.get_crate('HULLLIST')
    except IndexError as exc:
        utils.errlog("unable to read HULLLIST from " +
                     "hullfile {} - {}".format(hullfile, exc))
        out = "<!DOCTYPE html><html><head><title>INTERNAL ERROR</title>"
        out += "</head><body><p>Error reading HULLLIST block of "
        out += "hullfile={}\nreason=\n{}".format(hullfile, exc)
        out += "- see Doug!</p></body></html>"
        return 404, out

    for itervals in zip(cr.Master_Id.values,
                        cr.STATUS.values,
                        cr.BASE_STK.values,
                        cr.NVERTEX.values,
                        cr.EQPOS.values):

        midval, status, basestk, nvertex, eqpos = itervals
        assert midval in hull_store, \
            "Missing Master Id={}".format(midval)

        store = hull_store[midval]
        store['status'] = status
        store['basestack'] = basestk

        # The proposed hull
        #
        if nvertex > 0:
            wcs_orig = [eqpos_to_dict(eqpos, nvertex)]
        else:
            qafile = os.path.join(dname,
                                  'qa.{:03d}.v{}.fits'.format(midval,
                                                              revstr))
            if not os.path.isfile(qafile):
                utils.errlog("missing qafile {}".format(qafile))
                out = "<!DOCTYPE html><html><head><title>INTERNAL ERROR</title>"
                out += "</head><body><p>Missing qafile={}".format(qafile)
                out += "- see Doug!</p></body></html>"
                return 404, out

            try:
                qcr = pycrates.read_file(qafile + "[cols nvertex, eqpos]")
            except IOError as exc:
                utils.errlog("unable to read qafile {}".format(qafile))
                out = "<!DOCTYPE html><html><head><title>INTERNAL ERROR</title>"
                out += "</head><body><p>Unable to read "
                out += "qafile={}\nerror=\n{}".format(qafile, exc)
                out += "- see Doug!</p></body></html>"
                return 404, out

            wcs_orig = []
            for nvertex, eqpos in zip(qcr.NVERTEX.values,
                                      qcr.EQPOS.values):
                wcs_orig.append(eqpos_to_dict(eqpos, nvertex))

            qcr = None

        store['wcs_orig'] = wcs_orig

        # Has the user got their own version?
        polyname = utils.make_poly_name_json(ensemble, midval, revstr)
        polyfile = os.path.join(userdir, ensemble, polyname)
        if os.path.exists(polyfile):
            jcts = utils.read_json(polyfile)
            if jcts is None:
                utils.errlog("polyfile is unreadable: {}".format(polyfile))
        else:
            jcts = None

        if jcts is None:
            store['wcs'] = wcs_orig
        else:
            # Hack the dict to replace unicode keys as this breaks
            # things when converting back to JSON. Should this be
            # in read_json?
            #
            out = []
            for ps in jcts['polygons']:
                out.append({'ra': ps['ra'], 'dec': ps['dec']})

            store['wcs'] = out

    cr = None
    ds = None

    # What stacks are we interested in? That is, which stacks contain
    # hulls that contribute to this master hull?
    #
    cpts = hull_store[masterid]['components']
    stks = list(cpts.keys())

    # order by the stack number
    ordered_stks = sorted(stks, key=lambda stk: info['stackmap'][stk])

    enbands_cpt = {}
    for key, eband in ebands_by_component.items():
        # Had problems sending a pair via JSON, so convert to a string
        # for the key
        stk, cpt = key
        newkey = "{}.{}".format(stk, cpt)
        enbands_cpt[newkey] = eband

    # What about the stack-level hull polygons?
    #
    # Note that we include the center of the polygon,
    # calculated from its SKY coordinate values but converted
    # to RA and Dec. I know I have this code somewhere, but
    # "re-invent" it. Oh darn: at this point I don't easily
    # have the SKY coordinates, so going to hack it with the
    # equatorial coordinates (including a hack to deal with
    # ra=0/360).
    #
    stack_polys = {}
    for key in ebands_by_component.keys():

        stk, cpt = key
        regfile = os.path.join(dname,
                               'stack.{}.{}.v{}.reg'.format(stk,
                                                            cpt,
                                                            revstr))
        regstr = read_ds9_region(regfile)
        if regstr is None:
            utils.log("missing stack region file {}".format(regfile))
            out = "<!DOCTYPE html><html><head><title>INTERNAL ERROR</title>"
            out += "</head><body><p>Unable to read stack "
            out += "regfile={}".format(regfile)
            out += "- see Doug!</p></body></html>"
            return 404, out

        sra, sdec, slabel = regstr_to_coords(regstr)
        if slabel is None:
            # Use the stack number, not the full name, to save space
            # on the screen, if available.
            #
            try:
                stklbl = "{:03d}".format(stack_map[stk])
            except KeyError:
                # TODO: THIS SHOULD BE CONSIDERED AN ERROR
                stklbl = stk

            slabel = '{}.{:02d}'.format(stklbl, cpt)

        try:
            store = stack_polys[stk]
        except KeyError:
            store = []
            stack_polys[stk] = store

        # guestimates crossover limits

        sx = np.asarray(sra)
        sy = np.asarray(sdec)
        if (sx.min() < 20.0) and (sx.max() > 340):
            s0 = 400
        else:
            s0 = 0

        ra0, dec0 = utils.polygon_centroid(sx + s0, sy)
        ra0 -= s0

        shull = {'stack': stk, 'component': cpt,
                 'label': slabel,
                 'ra': sra, 'dec': sdec,
                 'ra0': ra0, 'dec0': dec0,
                 'regstr': regstr}
        try:
            shull['mancode'] = mancode_by_component[key]
        except KeyError:
            # this should not happen
            utils.log("key is missing in mancode_by_component: " +
                      "{}".format(key))
            # pass

        store.append(shull)

    # What are the point-source regions for the stacks?
    #
    psfs = {}
    for stk in stks:
        # hard code the name
        infile = os.path.join(xmdatdir, stk,
                              '{}N000_xmdat3.fits'.format(stk))

        # No VFS in the file name so can do an existence check
        if not os.path.exists(infile):
            continue

        try:
            cr = pycrates.read_file(infile)
        except IOError:
            utils.errlog("Unable to open xmdat file: {}".format(infile))
            continue

        # ra,dec in decmal degrees
        # r0,r1 in arcsec
        # ang in degrees
        #
        psfs[stk] = []
        for ra, dec, r0, r1, ang in zip(cr.ra.values,
                                        cr.dec.values,
                                        cr.psf_r0.values,
                                        cr.psf_r1.values,
                                        cr.psf_ang.values):
            psfs[stk].append({'ra': ra, 'dec': dec,
                              'r0': r0, 'r1': r1, 'angle': ang})

        cr = None

    # Ensure data needed by the templates is present (this should
    # be cleaned up)
    #
    for h in info['masters']:
        h['masterid_int'] = int(h['masterid'])

    detectors = ",".join(list(detectors))

    # We probably have this information, but convert it into a
    # form easily usable by the template.
    #
    # TODO: This overlaps with the data in the stack_polys dictionary,
    # which is keyed by the stack id.
    #
    components = []

    # Can not guarantee that the stacks are ordered lexicographically
    # within an ensemble, so sort them explicitly.
    #
    keylist = list(ebands_by_component.keys())
    keylist = sorted(keylist,
                     key=lambda key: (stack_map[key[0]], key[1]))

    for key in keylist:
        if mid_by_component[key] != masterid:
            continue

        name = "{:03d}.{:02d}".format(stack_map[key[0]], key[1])

        """
        TODO:
        these values need to be set up from the JSON file

        we need to have user/proposed values for the over-ridable values

        """

        # Ensure values are Python bools, not NumPy ones
        components.append({'name': name,
                           'masterid': int(masterid),
                           'eband': ebands_by_component[key],
                           'likelihood': likelihood_by_component[key],
                           'adjusted': mancode_by_component[key] > 0,
                           'svdqa':
                           bool(svdqa_by_component[key]),
                           'include_in_centroid':
                           bool(centroid_by_component[key]),
                           'mrg3rev': mrgrev_by_component[key]})

    # Do we know the actual size here?
    #
    if len(components) == 0:
        utils.log("filtered out all components!")
        out = "<!DOCTYPE html><html><head><title>INTERNAL ERROR</title>"
        out += "</head><body><p>Filtered out all components "
        out += "- see Doug!</p></body></html>"
        return 404, out

    # A lot of this data could be accessed via AJAX from the page
    # setup code, but leave as is for now.
    #
    # TODO: send in the "new" data
    return apply_template(env, 'masterhull.html',
                          {'ensemble': ensemble,
                           'revstr': revstr,
                           'mid': mid,
                           'masterid': masterid,
                           'npages': hull['npages'],
                           'ncpts': hull['ncpts'],
                           'detectors': detectors,
                           'components': components,
                           'info': info,
                           'hull': hull,
                           'is_latest': is_latest,
                           'ordered_stacks': ordered_stks,
                           'enbands_cpt': enbands_cpt,
                           'hull_store': hull_store,
                           'stack_polys': stack_polys,
                           'stack_psfs': psfs,
                           'username': os.getlogin()})


# Not sure what JS9 uses so cover some common types
#
mimetype_map = {
    '': 'text/plain',
    '.js': 'application/javascript',
    '.json': 'application/json',
    '.css': 'text/css',
    '.html': 'text/html',
    '.gif': 'image/gif',
    '.png': 'image/png',
    '.jpg': 'image/jpeg',
    '.jpeg': 'image/jpeg',
    '.wasm': 'application/wasm',
    # http://www.alienfactory.co.uk/articles/mime-types-for-web-fonts-in-bedsheet
    '.woff2': 'application/font-woff2',
    '.woff': 'application/font-woff',
    '.ttf': 'application/font-sfnt'
}


def get_mimetype(infile):
    """Return mime type for the file.

    This is purely based on the file name.

    Returns
    -------
    mimetype : str or None
        None is returned if the file is unknown.
    """

    suffix = os.path.splitext(infile)[1]
    try:
        return mimetype_map[suffix]
    except KeyError:
        return None


def send_file(obj, infile, mimetype, headers=None):
    """Send the contents of infile with the given mime type.

    Notes
    -----
    Can we stream this/does it make sense to do so?

    Should this just be a method of CHSHandler?

    Does it make sense to include Content-length?
    """

    try:
        cts = open(infile, 'rb').read()
    except IOError as exc:
        utils.log("error reading [{}]: {}".format(infile, exc))
        # could send something more useful
        obj.send_error(404)
        return

    obj.send_response(200)
    obj.send_header('Content-type', mimetype)
    obj.send_header('Content-length', len(cts))
    if headers is not None:
        for k, v in headers:
            obj.send_header(k, v)

    obj.end_headers()
    obj.write_bytes(cts)


def send_ensemble_hull_status(server,
                              datadir, userdir, ensemble, revision,
                              masterid):
    """Return - as JSON - the status information for this master hull.

    Parameters
    ----------
    server
        Used to send the response (404 or JSON).
    datadir : str
        The path to the directory containing the ensemble-level
        products.
    userdir : str
        The path to the directory containing the user's decisions.
    ensemble : str
        The ensemble (assumed to be valid at this point).
    revision : str
        This is assumed to be in 3-digit, 0 padded (left) form.
    masterid : int.
        The master id.

    Notes
    -----
    Should this include the component information? That is, the
    masterid and include_in_centroid values for each stack-level
    component in the master. At present this is being added in
    masterhull.js since we have the information there, but this is
    all rather muddled.
    """

    ensemble_status = utils.read_ensemble_status(datadir, userdir,
                                                 ensemble, revision)
    cts = utils.read_ensemble_hull_json(datadir, userdir,
                                        ensemble, masterid, revision)
    if cts is None:
        server.send_error(404)
        return

    cts['ensemble_status'] = ensemble_status

    # TODO: add in component information?

    server.send_as_json(cts)


class CHSHandler(BaseHTTPRequestHandler):

    def _set_headers(self, status=200):
        self.send_response(status)
        self.send_header('Content-type', 'text/html')
        self.end_headers()

    def write_bytes(self, contents, status=None):
        if status is not None:
            self._set_headers(status=status)
        self.wfile.write(contents)

    def write_contents(self, contents, status=None):
        # the encoding should be a no-op in Python 2.7
        self.write_bytes(contents.encode('utf-8'),
                         status=status)

    def do_GET(self):

        # This is going to be ugly routing code. Handle all
        # the special cases before dropping down to the
        # webassets directory (could special case these
        # but let's not for now).
        #
        # Note that the values are checked each call, rather than
        # cached. Let's see how this works.
        #
        context = self.server.context
        datadir = context['datadir']
        userdir = context['userdir']
        env = context['environment']
        upath = urlparse(self.path)
        path = upath.path

        # Strip the leading /
        if path != "":
            path = path[1:]

        if path in ["", "index.html"]:
            status, cts = create_index_page(env, datadir)
            self.write_contents(cts, status=status)
            return

        if path.startswith('evt3/'):
            self.get_evt3(path[5:])
            return

        if path.startswith('api/'):
            self.get_api(path[4:])
            return

        if path.startswith('img/'):
            self.get_image(path[4:])
            return

        # Is this an ensemble directory?
        #
        toks = path.split('/')
        if valid_ensemble(datadir, toks[0]):
            ntoks = len(toks)
            ensemble = toks[0]

            if ntoks == 1:
                status, cts = create_ensemble_index_page(env,
                                                         datadir,
                                                         userdir,
                                                         ensemble)
                self.write_contents(cts, status=status)
                return

            elif ntoks == 3:
                try:
                    revision = int(toks[1])
                    masterid = int(toks[2])
                except ValueError:
                    utils.errlog("Unable to parse revision/masterid " +
                                 "from [{}]".format(path))
                    self.send_error(404)
                    return

                # At the moment use the same data directory
                xmdatdir = context['xmdatdir']
                status, cts = create_master_hull_page(env,
                                                      datadir,
                                                      userdir,
                                                      datadir,
                                                      xmdatdir,
                                                      ensemble,
                                                      revision,
                                                      masterid)
                self.write_contents(cts, status=status)
                return

            else:
                utils.errlog("Invalid ensemble path " +
                             "[{}]".format(path))
                self.send_error(404)

        # Look in the webassets directory for this file
        #
        self.get_local_file(context['webdir'], path)

    def do_POST(self):

        upath = urlparse(self.path)
        path = upath.path

        # Strip the leading /
        if path != "":
            path = path[1:]

        try:
            savefunc = {'save/summary': utils.save_summary,
                        'save/ensemble': utils.save_ensemble,
                        'save/master': utils.save_master,
                        'save/masterpoly': utils.save_master_poly,
                        'save/component': utils.save_component,
                        'save/datatable': utils.save_datatable
                        }[path]
        except KeyError:
            utils.errlog("Unexpected POST path={}".format(path))
            self.send_error(404)
            return

        # for now require JSON
        ctype = self.headers['Content-Type'].split(';')
        if not ctype[0] == 'application/json':
            utils.errlog("Unexpected content-type: " +
                         "{}".format(ctype[0]))
            self.send_error(404)
            return

        clen = int(self.headers['Content-Length'])
        data = self.rfile.read(clen)
        try:
            jcts = json.loads(data)
        except Exception as exc:
            utils.errlog("Invalid JSON: {}".format(exc))
            self.send_error(404)
            return

        context = self.server.context
        userdir = context['userdir']

        try:
            savefunc(userdir, jcts)
        except Exception as exc:
            utils.errlog("Unable to save data: {}".format(exc))
            self.send_error(404)
            return

        self.write_contents("", status=200)

    def send_as_json(self, cts):
        """Return cts as JSON data to the caller."""

        try:
            rval = json.dumps(cts)
        except Exception as exc:
            utils.errlog("error converting JSON: {}".format(exc))
            self._set_headers(404)
            emsg = "<html><body><h1>ERROR</h1>" + \
                   "<p>Unable to convert data to JSON.</p>" + \
                   "<pre>{}</pre>".format(exc) + \
                   "</body></html>"
            self.write_contents(emsg, status=404)
            return

        self.send_response(200)
        self.send_header('Content-type', 'application/json')
        self.end_headers()
        self.write_contents(rval)

    def get_api(self, path):
        """Return the requested JSON data.

        path is assumed to have had the leading /api/ removed, so
        it does *not* begin with the '/' character.

        Supported are:
           summary
           ensxxxxxxx_xxx
           ensxxxxxxx_xxx/vvv/xxx

        """

        context = self.server.context
        datadir = context['datadir']
        userdir = context['userdir']
        if path == '':
            self.send_error(404)
            return

        elif path == 'summary':
            cts = get_data_summary(datadir, userdir)
            self.send_as_json(cts)
            return

        toks = path.split('/')
        ntoks = len(toks)

        ensemble = toks[0]
        if not valid_ensemble(datadir, ensemble):
            self.send_error(404)
            return

        if ntoks == 1:
            cts = utils.read_ensemble_json(datadir, userdir, ensemble)
            self.send_as_json(cts)
            return

        elif ntoks != 3:
            self.send_error(404)
            return

        # This is expected to already be in 001 form.
        revision = toks[1]

        try:
            masterid = int(toks[2])
        except ValueError:
            utils.errlog("Invalid master id: {}".format(toks[2]))
            self.send_error(404)
            return

        send_ensemble_hull_status(self, datadir, userdir, ensemble,
                                  revision, masterid)

    def get_image(self, path):
        """Return the image data.

        path is assumed to have had the leading /img/ removed, so
        it does *not* begin with the '/' character.

        Supported are:
           ensxxxxxxx_xxx/field.xxxxxxx_xxx.vxxx.png
           ensxxxxxxx_xxx/hull.xxxxxxx_xxx.xxx.pxxx.vxxx.<scale>.png

        Perhaps should not be so locked down.

        """

        if not path.endswith('.png'):
            self.send_error(404)
            return

        toks = path.split('/')
        ntoks = len(toks)
        if ntoks != 2:
            self.send_error(404)
            return

        if not toks[0].startswith('ens') or \
                not (toks[1].startswith('field.') or
                     toks[1].startswith('hull.')):
            self.send_error(404)
            return

        datadir = self.server.context['datadir']
        infile = os.path.join(datadir, '/'.join(toks))
        if not os.path.isfile(infile):
            self.send_error(404)
            return

        send_file(self, infile, 'image/png')

    def get_local_file(self, dirname, filename):
        """Return the file."""

        infile = os.path.join(dirname, filename)
        if not os.path.isfile(infile):
            self.send_error(404)
            return

        mimetype = get_mimetype(infile)
        if mimetype is None:
            utils.errlog("unknown JS9 suffix [{}]".format(infile))
            self.send_error(404)
            return

        send_file(self, infile, mimetype)

    def get_evt3(self, stack, use_gzip=False):
        """Return the evt3 file for the given stack.

        We handle both *.fits and *.fits.gz, so that
        large files can be uncompressed (to avoid running
        out of memory with the compression/decompression,
        mainly in other stages of the pipeline).

        use_gzip support is very-lightly tested; it is
        left in mainly for documentation.
        """

        # Need to get the version
        pat = os.path.join(self.server.context['evt3dir'],
                           "{}*fits*".format(stack))
        matches = glob.glob(pat)
        if len(matches) == 0:
            utils.errlog("failed glob={}".format(pat))
            self.send_error(404)
            return

        # Should pick the highest version; for now pick the first
        #
        if len(matches) > 1:
            utils.log("Multiple matches for {}".format(pat))

        infile = matches[0]

        if use_gzip:
            try:
                cts = open(infile, 'r').read()
            except IOError as exc:
                utils.errlog("error reading [{}]: ".format(infile) +
                             "{}".format(exc))

                # could send something more useful as a response?
                self.send_error(404)
                return

            self.send_response(200)
            self.send_header('Content-type', 'application/fits')
            self.send_header('Content-Encoding', 'gzip')
            # self.send_header('Content-length', len(cts))
            self.end_headers()

            with gzip.GzipFile(stack + ".fits",
                               compresslevel=9,
                               mode='wb',
                               fileobj=self.wfile) as fh:
                fh.write(cts)

        else:
            if infile.endswith('.gz'):
                hdrs = [('Content-Encoding', 'gzip')]
            else:
                hdrs = None

            send_file(self, infile, 'application/fits',
                      headers=hdrs)


def serve(userdir, webdir, datadir, templatedir,
          port=8070,
          xmdatdir='/data/L3/chs_master_match/input/xmdat3',
          evt3dir='/data/L3/chs_master_match/input/stkevt3'):

    for dirname in [userdir, webdir, datadir, templatedir,
                    evt3dir, xmdatdir]:
        if not os.path.isdir(dirname):
            raise IOError("Not a directory: {}".format(dirname))

    env = Environment(loader=FileSystemLoader(templatedir),
                      autoescape=select_autoescape(['html']))

    store = {'userdir': userdir,
             'webdir': webdir,
             'datadir': datadir,
             'evt3dir': evt3dir,
             'xmdatdir': xmdatdir,
             'environment': env}

    server_address = ('', port)
    httpd = CHSServer(store, server_address, CHSHandler)
    utils.log("Web assets: {}".format(webdir))
    utils.log("Starting server on http://localhost:{}/".format(port))
    httpd.serve_forever()


help_str = "Run the CHS QA GUI (web pages)."


if __name__ == "__main__":

    import argparse

    parser = argparse.ArgumentParser(description=help_str,
                                     prog=sys.argv[0])

    parser.add_argument("datadir",
                        help="The directory containing the mhull and image files")
    parser.add_argument("--port", type=int, default=8070,
                        help="The port number to run the QA server on: default %(default)i")

    parser.add_argument("--xmdatdir",
                        default="/data/L3/chs_master_match/input/xmdat3",
                        help="The xmdat3 directory: default %(default)s")
    parser.add_argument("--evt3dir",
                        default="/data/L3/chs_master_match/input/stkevt3",
                        help="The stkevt3 directory: default %(default)s")

    args = parser.parse_args(sys.argv[1:])

    datadir = os.path.normpath(os.path.abspath(args.datadir))

    # The assets for the web server are obtained from the location
    # of the script (../webassets/).
    #
    thisdir = os.path.dirname(os.path.realpath(__file__))
    webassets = os.path.join(thisdir, '../webassets')
    webassets = os.path.normpath(webassets)

    # Want a directory that can be used in testing.
    #
    templatedir = os.path.join(webassets, 'templates')

    # Store the user's files in the current working directory
    #
    userdir = os.path.normpath(os.getcwd())

    # Validation check: ensure that the user directory is *not*
    # within the data directory, to avoid over-writing files.
    #
    cdir = os.path.commonprefix([datadir, userdir])

    def stripify(d):
        if d.endswith('/'):
            return d[:-1]
        else:
            return d

    if stripify(cdir) == stripify(datadir) or \
            stripify(datadir) == stripify(userdir):
        raise IOError("Userdir {} is a ".format(userdir) +
                      "sub-directory of " +
                      "Datadir {}".format(datadir))

    serve(userdir=userdir, webdir=webassets,
          datadir=datadir, templatedir=templatedir,
          port=args.port,
          xmdatdir=args.xmdatdir,
          evt3dir=args.evt3dir)
