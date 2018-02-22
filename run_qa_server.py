#!/usr/bin/env python

"""

Usage:

  ./run_qa_server.py datadir [port]

The default port is 8070. The datadir is the location of the
master-hull data (the location of the ens<> directories).

The webserver assets are found relative to the location of the
script.

Data will be written to the working directory (or it will eventually
do so).

Things to do

 - hull review page: regions could be added to JS when the page is created, which would
   avoid some of the callbacks needed
 - should draw all stack hulls in JS9 window

  These two are almost done; need to remove or change the behavior of the
  reload region button, then can remove/clean up the region endpoint
  and need to send in the hull centeres for the master hulls

  Need to start labelling the master hulls (for when there are qa cases)
  so that can delete a single hull.

  How to delete the convex hull after adjustment? Maybe don't do this
  ambiently, but wait for the user to save the hull, in which case the
  display can be updated at that time (or perform some other action).

 - when the master hull is changed, should change the JS data representation and then use
   that to update other JS9 displays.

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

It is Python 2.7 only at present.

See also https://pymotw.com/2/BaseHTTPServer/

"""

import glob
import gzip
import json
import os
import sys

from six.moves.BaseHTTPServer import BaseHTTPRequestHandler, \
    HTTPServer

import urlparse

import numpy as np

from jinja2 import Environment, FileSystemLoader, select_autoescape, \
    TemplateNotFound

try:
    import pycrates
except ImportError:
    sys.stderr.write("Needs crates!\n")
    sys.exit(1)


def log(msg, level='LOG'):
    """Log a message.

    Parameters
    ----------
    msg: str
        The message
    """

    print("{}: {}".format(level, msg))


def dbg(msg):
    """Log a debug message.

    Parameters
    ----------
    msg: str
        The message
    """

    log(msg, level='DEBUG')


def warn(msg):
    """Log a warning message.

    Parameters
    ----------
    msg: str
        The message
    """

    log(msg, level='WARNING')


def errlog(msg):
    """Log an error message.

    Parameters
    ----------
    msg: str
        The message
    """

    log(msg, level='ERROR')


class CHSServer(HTTPServer):
    """Allow information to be passed to the handler.

    The store argument should be a dictionary with
    keywords:
       datadir  - the directory containing the ensemble/master hull
                  information (this is created by the CHS
                  master-match pipeline)
       evt3dir  - the location of the stack event files
       userdir  - the location for storing user information
       webdir   - the location of the web assets
       environment - the Jinja2 environment

    These must all be existing directories.
    """

    def __init__(self, store, *args, **kwargs):
        context = {}
        for k in ['datadir', 'evt3dir', 'userdir', 'webdir']:
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

    try:
        int(ensemble[3:-4])
        int(ensemble[-3:])
    except ValueError:
        return False

    path = os.path.join(datadir, ensemble)
    return os.path.isdir(path)


def read_json(infile):
    """Return data stored as JSON in the input file.

    Parameters
    ----------
    infile : str
        The full path to the JSON file.

    Returns
    -------
    ans : dict or array or None
        The contents, or None if there was an error.
    """

    try:
        with open(infile, 'r') as fh:
            cts = fh.read()
            jcts = json.loads(cts)

    except Exception as exc:
        errlog("error reading JSON from {}\n{}".format(infile, exc))
        return None

    return jcts


def read_ensemble_hull_json(datadir, ensemble, mid, revision):
    """Return JSON-stored data for the masters.

    Parameters
    ----------
    datadir : str
        The path to the directory containing the ensemble-level
        products.
    ensemble : str
        The ensemble.
    mid : int
        The master id.
    revision : str
        Expected to be 001, ...

    Returns
    -------
    ans : dict or None
        The contents, or None if there was an error.

    """

    infile = os.path.join(datadir, ensemble,
                          'hull.{}.{:03d}.v{}.json'.format(ensemble,
                                                           mid,
                                                           revision))
    jcts = read_json(infile)
    if jcts is None:
        return None

    # raise NotImplementedError()
    return jcts


def read_ensemble_json(datadir, ensemble):
    """Return JSON data from the ensemble.

    The returned structure contains all the versions for this
    ensemble. The "current" version is taken to be the highest
    version number.

    Parameters
    ----------
    datadir : str
        The path to the directory containing the ensemble-level
        products.
    ensemble : str
        The ensemble.

    Returns
    -------
    ans : dict or None
        The contents, or None if there was an error.

    Notes
    -----

    state['latest_version'] => str (001, 002, ..)
    state['versions'] => dict, keys are version
          ['nstacks']
          ['nmasters']  == len(masters)
          ['masters'] = list of dict
                  ['masterid'] => str (001, ...)
                  ['ncpts'] == len(cpts)
                  ---- ['cpts'] => list of ?
                  ['useraction'] => str, can be ''
                  ['lastmodified'] =>

hmmm, the JSON data has {"status": "todo", "stackmap": {"acisfJ1705367m403832_001": 3, "acisfJ1704041m414416_001": 1, "acisfJ1702545m412821_001": 0, "acisfJ1705559m410515_001": 4, "acisfJ1704448m410953_001": 2}, "nmasters": 1, "name": "ens0000900_001", "lastmodified": "", "ncpts": 2, "usernotes": "", "nstacks": 2, "revision": "001"}


    """

    pat = "field.{}.*.json".format(ensemble)
    inpat = os.path.join(datadir, ensemble, pat)
    matches = glob.glob(inpat)
    if len(matches) == 0:
        errlog("no field.json files found for ensemble " +
               "{} - {}".format(ensemble, inpat))
        return None

    store = {'versions': {}}
    for match in matches:
        jcts = read_json(match)
        if jcts is None:
            continue

        try:
            v = jcts['revision']
        except KeyError:
            log("missing revision keyword in {}".format(match))
            continue

        # this should not happen, so do not worry too much about the
        # error handler
        if v in store['versions']:
            log("multiple revision={} in {}".format(v, match))
            continue

        hulls = []
        for mid in range(1, jcts['nmasters'] + 1):
            hull = read_ensemble_hull_json(datadir, ensemble, mid, v)
            if hull is None:
                # if there'sa problem reading in a single hull, then
                # bail out for the whole thing
                return None

            hulls.append(hull)

        if len(hulls) == 0:
            log("no masters for revision {} in {}".format(v, match))
            continue

        if 'masters' in jcts:
            log("overwriting masters setting in {}".format(match))

        jcts['masters'] = hulls
        store['versions'][v] = jcts

    revs = list(store['versions'].keys())
    if len(revs) == 0:
        errlog("no JSON data read from ensemble: " +
               "{} {}".format(datadir, ensemble))
        return None

    revs = sorted(revs, key=int, reverse=True)
    store['latest_version'] = revs[0]

    return store


def parse_datadir(datadir):
    """Extract useful information from the diretory.

    Parameters
    ----------
    datadir : str
        The path to the directory containing the ensemble-level
        products.

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
        errlog("no ensemble dirs found in {}".format(datadir))
        return None

    # HACK: use a subset of ensembles to save time
    #
    """
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
        jcts = read_ensemble_json(ensembledir + "/../", ensemble)
        if jcts is None:
            continue

        v = jcts['latest_version']
        info = jcts['versions'][v]
        e = info['name']
        store[e] = info

    # return in ensemble order
    keys = sorted(store.keys())
    return [store[k] for k in keys]


def parse_ensembledir(datadir, ensemble):
    """Extract useful information from the diretory.

    This extends read_ensemble_json to add in the per-master
    information.

    Parameters
    ----------
    datadir : str
        The path to the directory containing the ensemble-level
        products.
    ensemble : str
        The ensemble.

    Returns
    -------
    state : dict or None
        A dictionary or None, if no data can be found there.

    Notes
    -----
    At present the contents of datadir are highly structured;
    let's see if this will need relaxing.
    """

    jcts = read_ensemble_json(datadir, ensemble)
    if jcts is None:
        return None

    # TODO: STUFF
    return jcts


def get_data_summary(datadir):
    """Return the data about all the available ensembles.

    Parameters
    ----------
    datadir : str
        The path to the directory containing the ensemble-level
        products.

    Returns
    -------
    summary : dict or None
    """

    state = parse_datadir(datadir)
    if state is None:
        return state

    # split up the ensembles
    #
    # QUESTION: what are the valid status values?
    #
    out = {'todos': [],
           'reviews': [],
           'completed': []}
    for ens in state:
        if ens['status'] == 'completed':
            key = 'completed'
        elif ens['status'] == 'review':
            key = 'reviews'
        else:
            key = 'todos'

        out[key].append(ens)

    return out


def get_data_ensemble(datadir, ensemble):
    """Return the data about this ensemble.

    Parameters
    ----------
    datadir : str
        The path to the directory containing the ensemble-level
        products.
    ensemble : str
        The ensemble id.

    Returns
    -------
    summary : dict or None
    """

    state = parse_ensembledir(datadir, ensemble)
    return state


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
                        errlog("expected coordsys, found " +
                               "[{}] in {}".format(l, infile))
                        return None

                    csys = l
                    continue

                if not l.startswith('polygon('):
                    errlog("expected polygon(...), found " +
                           "[{}] in {}".format(l, infile))
                    return None

                return "{}; {}".format(csys, l)

    except IOError as exc:
        errlog("can not read region file:" +
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


def find_master_center(infile, masterid):
    """Return the approximate center of the master hull.

    Parameters
    ----------
    infile : str
        The name of the file containing the master hull.
    masterid : int
        The matching Master_Id value.

    Returns
    -------
    ans : None or (ra, dec)
        The approximate center of the hull. There is no accounting
        for spherical projection, and isn't the polygon centroid.
        If the hull crosses ra=0/360 then the answer will be wrong.
        Returns None if this is a QA case or there is a problem
        reading the file.

    """

    fname = "{}[SRCLIST][Master_Id={}]".format(infile, masterid)
    cr = pycrates.read_file(fname)
    nrows = cr.get_nrows()
    if nrows == 0:
        errlog("no row matching {}".format(fname))
        return None

    elif nrows > 1:
        errlog("multiple rows matching {}".format(fname))
        return None

    # NOTE: QA cases have NVERTEX=0
    nv = cr.NVERTEX.values[0]
    if nv < 3:
        errlog("nvertex={} for {}".format(nv, fname))
        return None

    eqpos = cr.EQPOS.values[0, :, :nv]
    ra_mid = (eqpos[0].min() + eqpos[0].max()) / 2.0
    dec_mid = (eqpos[1].min() + eqpos[1].max()) / 2.0
    return ra_mid, dec_mid


def get_data_master(datadir, rawdir, ensemble, masterid):
    """Return the data about this master hull.

    Parameters
    ----------
    datadir : str
        The path to the directory containing the ensemble-level
        products.
    rawdir : str
        The path to the actual master hull data (fits and region files)
    ensemble : str
        The ensemble id.
    masterid : int

    Returns
    -------
    summary : dict or None
    """

    # Need the latest version of the master hull. Now
    # read_ensemble_json does too much, but let's not worry about
    # that here.
    state = read_ensemble_json(datadir, ensemble)
    if state is None:
        return None

    try:
        version = state['latest_version']
        masters = state['versions'][version]['masters']
    except KeyError as exc:
        errlog("ensemble json missing key {}".format(exc))
        return None

    mid = "{:03d}".format(masterid)
    hull = None
    for m in masters:
        if m['masterid'] == mid:
            hull = m
            break

    if hull is None:
        log("missing hull {}".format(mid))
        return None

    # At this point we are happy that ensemble/masterid point to a
    # "valid" hull. Not sure the above was strictly necessary
    # as going to the following.
    #
    # Now need to parse the actual hull data: at the moment use
    # the actual data but maybe this should be in the JSON?
    #
    dname = os.path.join(rawdir, ensemble)
    if not os.path.isdir(dname):
        errlog("missing rawdir {}".format(dname))
        return None

    hullfile = os.path.join(dname,
                            'master_hulls.' +
                            '{}.v{}.fits'.format(ensemble, version))
    if not os.path.isfile(hullfile):
        errlog("missing hullfile {}".format(hullfile))
        return None

    infile = "{}[master_id={}][cols stackid, component]".format(hullfile,
                                                                masterid)
    cr = pycrates.read_file(infile)
    nrows = cr.get_nrows()
    if nrows == 0:
        errlog("no hull data from {}".format(infile))
        return None

    # Just warn, nbut do not error out, if there's a difference
    #
    if nrows != hull['ncpts']:
        warn("expected {} but got {} from {}".format(hull['ncpts'],
                                                     nrows,
                                                     infile))

    stackmap = {}
    nstacks = cr.get_key_value('STKIDNUM')
    if nstacks is not None:
        for i in range(nstacks):
            k = '{:03d}'.format(i)
            stackmap[cr.get_key_value('STKID{}'.format(k))] = k

    stacks = cr.STACKID.values.copy()
    components = cr.COMPONENT.values.copy()
    cr = None

    infile = os.path.join(dname,
                          'master.{}.v{}.reg'.format(masterid,
                                                     version))

    out = {}

    # The conversion to coordinates from the regstr is ugly
    # but an easy hack given the current code and data setup
    # (e.g. we don't guarantee the stack-level mrgsrc file
    # needed to get the coordinates of the stack-level hulls).
    #

    # Allow the master region to be undefined (e.g. QA case)
    regstr = read_ds9_region(infile)
    if regstr is not None:
        rra, rdec, rlabel = regstr_to_coords(regstr)
        if rlabel is None:
            rlabel = mid

        out['master'] = {'label': rlabel,
                         'ra': rra,
                         'dec': rdec,
                         'regstr': regstr}

    # For now have the center as a separate bit of logic to the
    # region, even though we should have both or neither, not
    # either.
    #
    infile = os.path.join(dname,
                          'master_hulls.{}.v{}.fits'.format(ensemble,
                                                            version))
    hull_cen = find_master_center(infile, int(masterid))
    if hull_cen is not None:
        out['center'] = {'ra': hull_cen[0], 'dec': hull_cen[1]}

    out['stacks'] = []
    for stk, cpt in zip(stacks, components):
        infile = os.path.join(dname,
                              'stack.{}.{}.v{}.reg'.format(stk,
                                                           cpt,
                                                           version))
        regstr = read_ds9_region(infile)
        if regstr is None:
            continue

        sra, sdec, slabel = regstr_to_coords(regstr)
        if slabel is None:
            # Use the stack number, not the full name, to save space
            # on the screen, if available.
            #
            try:
                stklbl = stackmap[stk]
            except KeyError:
                stklbl = stk

            slabel = '{}.{}'.format(stklbl, cpt)

        out['stacks'].append({'stack': stk, 'component': cpt,
                              'label': slabel,
                              'ra': sra,
                              'dec': sdec,
                              'regstr': regstr})
    return out


def get_master_hull_regions(ensemble, masterid, datadir, rawdir):
    """Return the master and stack hulls.

    Parameters
    ----------
    ensemble : str
    masterid : int
    datadir : str
    rawdir : str

    Returns
    -------
    ans : dict or None

    """

    return get_data_master(datadir, rawdir, ensemble, masterid)


def apply_template(env, tmplname, args):

    try:
        tmpl = env.get_template(tmplname)
    except TemplateNotFound as exc:
        errlog("No template for " + tmplname)
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


def create_ensemble_index_page(env, datadir, ensemble):
    """Create the top-level ensemble page.

    Parameters
    ----------
    datadir : str
        The path to the directory containing the ensemble-level
        products.
    ensemble : str
        The ensemble id.

    Returns
    -------
    status, page : int, str
        The response status and HTML contents
    """

    state = parse_ensembledir(datadir, ensemble)
    if state is None:
        out = "<!DOCTYPE html><html><head><title>ERROR</title>"
        out += "</head><body><p>There was an error when processing "
        out += "the directory {} for ".format(datadir)
        out += "ensemble {}".format(ensemble)
        out += "</body></html>"
        return 404, out

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
                            datadir, rawdir,
                            ensemble, revision, masterid):
    """Create the review page for a master hull.

    This is not ideal, since the master id depends on the version;
    i.e. v001 may have 2 masters but v002 may have 1 or 3. In fact,
    the master id number is not very stable - it can end up referring
    to different master hulls in different versions.

    Parameters
    ----------
    datadir : str
        The path to the directory containing the ensemble-level
        products.
    rawdir : str
        The path to the actual master hull data (fits and region files)
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

    mid = "{:03d}".format(masterid)
    revstr = "{:03d}".format(revision)

    state = parse_ensembledir(datadir, ensemble)
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
        log("hulls = {}".format(hulls))
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
        errlog("missing rawdir {}".format(dname))
        out = "<!DOCTYPE html><html><head><title>INTERNAL ERROR</title>"
        out += "</head><body><p>Missing dir={}".format(dname)
        out += + "- see Doug!</p></body></html>"
        return 404, out

    hullfile = os.path.join(dname,
                            'master_hulls.' +
                            '{}.v{}.fits'.format(ensemble, revstr))
    if not os.path.isfile(hullfile):
        errlog("missing hullfile {}".format(hullfile))
        out = "<!DOCTYPE html><html><head><title>INTERNAL ERROR</title>"
        out += "</head><body><p>Missing hullfile={}".format(hullfile)
        out += + "- see Doug!</p></body></html>"
        return 404, out

    # From the hull file,
    #
    # *) SRCMATCH block
    #
    # find those stacks associated with a master source, and the
    # energy band.
    #
    # *) SRCLIST block
    #
    # find the state of each master hull and, where possible,
    # read in the vertices.
    #
    try:
        ds = pycrates.CrateDataset(hullfile, mode='r')
    except IOError as exc:
        errlog("unable to read hullfile {} - {}".format(hullfile, exc))
        out = "<!DOCTYPE html><html><head><title>INTERNAL ERROR</title>"
        out += "</head><body><p>Error reading "
        out += "hullfile={}\nreason=\n{}".format(hullfile, exc)
        out += + "- see Doug!</p></body></html>"
        return 404, out

    try:
        cr = ds.get_crate('SRCMATCH')
    except IndexError as exc:
        errlog("unable to read SRCMATCH from hullfile {} - {}".format(hullfile, exc))
        out = "<!DOCTYPE html><html><head><title>INTERNAL ERROR</title>"
        out += "</head><body><p>Error reading SRCMATCH block of "
        out += "hullfile={}\nreason=\n{}".format(hullfile, exc)
        out += + "- see Doug!</p></body></html>"
        return 404, out

    # We have aleady got this information, but recreate it
    #
    nstks = cr.get_key_value('STKIDNUM')
    if nstks is None:
        errlog("missing STKIDNUM keyword in hullfile {} - {}".format(hullfile))
        out = "<!DOCTYPE html><html><head><title>INTERNAL ERROR</title>"
        out += "</head><body><p>No STKIDNUM keyword in "
        out += "hullfile={}\n".format(hullfile)
        out += + "- see Doug!</p></body></html>"
        return 404, out

    stack_map = {}
    for i in range(nstks):
        key = "STKID{:03d}".format(i)
        val = cr.get_key_value(key)
        if val is None:
            errlog("missing {} keyword in hullfile {} - {}".format(key, hullfile))
            out = "<!DOCTYPE html><html><head><title>INTERNAL ERROR</title>"
            out += "</head><body><p>No {} keyword in ".format(key)
            out += "hullfile={}\n".format(hullfile)
            out += + "- see Doug!</p></body></html>"
            return 404, out

        stack_map[val] = i

    # Do not need the eband info to be indexed by master id
    hull_store = {}
    ebands_by_component = {}
    mancode_by_component = {}
    for vals in zip(cr.Master_Id.values,
                    cr.STACKID.values,
                    cr.COMPONENT.values,
                    cr.EBAND.values,
                    cr.MAN_CODE.values):

        midval, stackid, cpt, eband, mancode = vals

        # NOTE: mancode is encoded "strangely", so remove the
        #       extra dimension
        if mancode.shape == (1,):
            mancode = mancode[0]

        # mancode is also stored as np.uint8 which seems to be problematic
        # to serialize to JSON, so explicitly convert to a "normal" int
        #
        mancode = int(mancode)

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
        assert key not in ebands_by_component, "key={}".format(key)
        assert key not in mancode_by_component, "key={}".format(key)
        ebands_by_component[key] = eband
        mancode_by_component[key] = mancode

    cr = None

    try:
        cr = ds.get_crate('SRCLIST')
    except IndexError as exc:
        errlog("unable to read SRCLIST from hullfile {} - {}".format(hullfile, exc))
        out = "<!DOCTYPE html><html><head><title>INTERNAL ERROR</title>"
        out += "</head><body><p>Error reading SRCLIST block of "
        out += "hullfile={}\nreason=\n{}".format(hullfile, exc)
        out += + "- see Doug!</p></body></html>"
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

        if nvertex > 0:
            wcs = [eqpos_to_dict(eqpos, nvertex)]
        else:
            qafile = os.path.join(dname,
                                  'qa.{}.v{}.fits'.format(mid,
                                                          revstr))
            if not os.path.isfile(qafile):
                errlog("missing qafile {}".format(qafile))
                out = "<!DOCTYPE html><html><head><title>INTERNAL ERROR</title>"
                out += "</head><body><p>Missing qafile={}".format(qafile)
                out += "- see Doug!</p></body></html>"
                return 404, out

            try:
                qcr = pycrates.read_file(qafile + "[cols nvertex, eqpos]")
            except IOError as exc:
                errlog("unable to read qafile {}".format(qafile))
                out = "<!DOCTYPE html><html><head><title>INTERNAL ERROR</title>"
                out += "</head><body><p>Unable to read "
                out += "qafile={}\nerror=\n{}".format(qafile, exc)
                out += "- see Doug!</p></body></html>"
                return 404, out

            wcs = []
            for nvertex, eqpos in zip(qcr.NVERTEX.values,
                                      qcr.EQPOS.values):
                wcs.append(eqpos_to_dict(eqpos, nvertex))

            qcr = None

        store['wcs'] = wcs

    cr = None
    ds = None

    # What stacks are we interested in? That is, which stacks contain
    # hulls that contribute to this master hull?
    #
    cpts = hull_store[masterid]['components']
    stks = list(cpts.keys())

    # order by the stack number
    ordered_stks = sorted(stks, key=lambda stk: info['stackmap'][stk])

    # Restrict the band information just to those stack-level hulls
    # that form the master hull. This can mean that the display
    # is not optimised for stack-level hulls from another master
    # source, bit these are informational only.
    #
    # The stack band information is only used for ACIS stacks,
    # and is a combination of the bands for the individual
    # components (at least at present).
    #
    minband = {'b': 500, 'u': 200, 's': 500, 'm': 1200, 'h': 2000}
    maxband = {'b': 7000, 'u': 500, 's': 1200, 'm': 2000, 'h': 7000}

    stack_band = {}
    for stk in stks:
        if stk.startswith('h'):
            continue

        elo = 10000
        ehi = 0
        for cpt in cpts[stk]:
            eband = ebands_by_component[(stk, cpt)]
            elo = min(elo, minband[eband])
            ehi = max(ehi, maxband[eband])

        stack_band[stk] = 'energy >= {} && energy < {}'.format(elo,
                                                               ehi)

    # What about the stack-level hull polygons?
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
            log("missing stack region file {}".format(regfile))
            out = "<!DOCTYPE html><html><head><title>INTERNAL ERROR</title>"
            out += "</head><body><p>Unable to read stack "
            out += "regfile={}".format(regfile)
            out += + "- see Doug!</p></body></html>"
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

        shull = {'stack': stk, 'component': cpt,
                 'label': slabel,
                 'ra': sra, 'dec': sdec,
                 'regstr': regstr}
        try:
            shull['mancode'] = mancode_by_component[key]
        except KeyError:
            # this should not happen
            log("key is missing in mancode_by_component: {}".format(key))
            # pass

        store.append(shull)

    # Ensure data needed by the templates is present (this should
    # be cleaned up)
    #
    for h in info['masters']:
        h['masterid_int'] = int(h['masterid'])

    return apply_template(env, 'masterhull.html',
                          {'ensemble': ensemble,
                           'revstr': revstr,
                           'mid': mid,
                           'masterid': masterid,
                           'npages': hull['npages'],
                           'ncpts': hull['ncpts'],
                           'info': info,
                           'hull': hull,
                           'is_latest': is_latest,
                           'ordered_stacks': ordered_stks,
                           # 'stack_band': json.dumps(stack_band),
                           # 'hull_store': json.dumps(hull_store),
                           # 'stack_polys': json.dumps(stack_polys),
                           'stack_band': stack_band,
                           'hull_store': hull_store,
                           'stack_polys': stack_polys,
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
        cts = open(infile, 'r').read()
    except IOError as exc:
        log("error reading [{}]: {}".format(infile, exc))
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
    obj.wfile.write(cts)


class CHSHandler(BaseHTTPRequestHandler):

    def _set_headers(self, status=200):
        self.send_response(status)
        self.send_header('Content-type', 'text/html')
        self.end_headers()

    def write_contents(self, status, contents):
        self._set_headers(status=status)
        self.wfile.write(contents)

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
        env = context['environment']
        upath = urlparse.urlparse(self.path)
        path = upath.path

        # Strip the leading /
        if path != "":
            path = path[1:]

        if path in ["", "index.html"]:
            status, cts = create_index_page(env, datadir)
            self.write_contents(status, cts)
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

        if path.startswith('regions/'):
            self.get_regions(path[8:])
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
                                                         ensemble)
                self.write_contents(status, cts)
                return

            elif ntoks == 3:
                try:
                    revision = int(toks[1])
                    masterid = int(toks[2])
                except ValueError:
                    errlog("Unable to parse revision/masterid " +
                           "from [{}]".format(path))
                    self.send_error(404)
                    return

                # At the moment use the same data directory
                status, cts = create_master_hull_page(env,
                                                      datadir,
                                                      datadir,
                                                      ensemble,
                                                      revision,
                                                      masterid)
                self.write_contents(status, cts)
                return

            else:
                errlog("Invalid ensemble path [{}]".format(path))
                self.send_error(404)

        # Look in the webassets directory for this file
        #
        self.get_local_file(context['webdir'], path)

    def do_POST(self):
        clen = int(self.headers['Content-Length'])
        data = self.rfile.read(clen)
        self._set_headers()
        self.wfile.write("<html><body><h1>POST!</h1><pre>")
        self.wfile.write(data)
        self.wfile.write("</pre></body></html>")

    def send_as_json(self, cts):
        """Return cts as JSON data to the caller."""

        try:
            rval = json.dumps(cts)
        except Exception as exc:
            errlog("error converting JSON: {}".format(exc))
            self._set_headers(404)
            self.wfile.write("<html><body><h1>ERROR</h1>")
            self.wfile.write("<p>Unable to convert data to JSON.</p>")
            self.wfile.write("<pre>{}</pre>".format(exc))
            self.wfile.write("</body></html>")
            return

        self.send_response(200)
        self.send_header('Content-type', 'application/json')
        self.end_headers()
        self.wfile.write(rval)

    def get_api(self, path):
        """Return the requested JSON data.

        path is assumed to have had the leading /api/ removed, so
        it does *not* begin with the '/' character.

        Supported are:
           summary
           ensxxxxxxx_xxx
           ensxxxxxxx_xxx/xxx   NOT USED

        Note that, as coded, the ensxxx/yyy version can not work,
        so need to see if actually needed.

        """

        context = self.server.context
        datadir = context['datadir']
        if path == '':
            self.send_error(404)
            return

        elif path == 'summary':
            cts = get_data_summary(datadir)
            self.send_as_json(cts)
            return

        toks = path.split('/')
        ntoks = len(toks)
        if ntoks > 2:
            self.send_error(404)
            return

        ensemble = toks[0]
        if not valid_ensemble(datadir, ensemble):
            self.send_error(404)
            return

        if ntoks == 1:
            cts = get_data_ensemble(datadir, ensemble)
            self.send_as_json(cts)
            return

        try:
            masterid = int(toks[1])
        except ValueError:
            self.send_error(404)

        # I DO NOT BELIEVE THIS PART OF THE CODE IS REACHED
        errlog("UNEXPECTED CODE EXECUTION")
        cts = get_data_master(datadir, ensemble, masterid)
        self.send_as_json(cts)

    def get_regions(self, path):
        """Return the requested region data as JSON.

        Supported are:
           ensemble/ensxxxxxxx_xxx/master_hull_id

        where master_hull_id may be 001 or 1

        """

        toks = path.split('/')
        if len(toks) != 3 or toks[0] != 'ensemble':
            self.send_error(404)
            return

        ensemble = toks[1]
        try:
            masterid = int(toks[2])
        except ValueError:
            log("invalid master id: [{}]".format(toks[2]))
            self.send_error(404)

        datadir = self.server.context['datadir']
        if not valid_ensemble(datadir, ensemble):
            self.send_error(404)
            return

        regs = get_master_hull_regions(ensemble, masterid,
                                       datadir, datadir)
        if regs is None:
            self.send_error(404)
            return None

        self.send_as_json(regs)

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
            errlog("unknown JS9 suffix [{}]".format(infile))
            self.send_error(404)
            return

        send_file(self, infile, mimetype)

    def get_evt3(self, stack, use_gzip=False):
        """Return the evt3 file for the given stack.

        Let's see how gzip-encoding works: if the input
        is gzip-encoded then there's less work for the server!
        """

        # Need to get the version
        pat = os.path.join(self.server.context['evt3dir'],
                           "{}*fits*".format(stack))
        matches = glob.glob(pat)
        if len(matches) == 0:
            errlog("failed glob={}".format(pat))
            self.send_error(404)
            return

        # Should pick the highest version; for now pick the first
        #
        infile = matches[0]

        if use_gzip:
            try:
                cts = open(infile, 'r').read()
            except IOError as exc:
                errlog("error reading [{}]: ".format(infile) +
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
          evt3dir='/data/L3/chs_master_match/input/stkevt3'):

    for dirname in [userdir, webdir, datadir, evt3dir, templatedir]:
        if not os.path.isdir(dirname):
            raise IOError("Not a directory: {}".format(dirname))

    env = Environment(loader=FileSystemLoader(templatedir),
                      autoescape=select_autoescape(['html']))

    store = {'userdir': userdir,
             'webdir': webdir,
             'datadir': datadir,
             'evt3dir': evt3dir,
             'environment': env}

    server_address = ('', port)
    httpd = CHSServer(store, server_address, CHSHandler)
    log("Starting server on http://localhost:{}/".format(port))
    httpd.serve_forever()


def usage(progName):
    sys.stderr.write("Usage: {} datadir [port]\n".format(sys.argv[0]))
    sys.stderr.write("\nport should be an integer and defaults to ")
    sys.stderr.write("8070 if not given.\n")
    sys.exit(1)

if __name__ == "__main__":

    nargs = len(sys.argv)
    if nargs == 2:
        port = 8070
    elif nargs == 3:
        try:
            port = int(sys.argv[2])
        except ValueError:
            usage(sys.argv[0])

    else:
        usage(sys.argv[0])

    datadir = os.path.abspath(sys.argv[1])

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
    userdir = os.getcwd()

    serve(userdir=userdir, webdir=webassets,
          datadir=datadir, templatedir=templatedir,
          port=port)
