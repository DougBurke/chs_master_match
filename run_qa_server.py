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

It is Python 2.7 only at present.

See also https://pymotw.com/2/BaseHTTPServer/

"""

import glob
import gzip
import json
import os
import sys
import time

import six
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

import chs_utils as utils


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


def setup_user_setting(store, key):
    """Add in user-level version.

    This is highly-specialized. It changes the key value to
    be a dictionary with 'proposed' and 'user' keywords,
    storing the current value under 'proposed' and setting
    'user' to None.

    Parameters
    ----------
    store : dict
        A dictiionary which is assumed to contain the supplied
        key, and the value of the key is a string.
    key : dict_key
        The key value.

    Returns
    -------
    flag : bool
        If False then the stored value was nto a string so nothing
        has been changed. It is expected that this will cause
        down stream to trigger an error handler.

    """

    try:
        v = store[key]
    except KeyError:
        warn("NO {} field".format(key))
        v = ''

    if not isinstance(v, six.string_types):
        errlog("{} is not a string but {}".format(key, v))
        return False

    store[key] = {'proposed': v, 'user': None}
    return True


def read_ensemble_hull_json(datadir, userdir,
                            ensemble, mid, revision):
    """Return JSON-stored data for the masters.

    Parameters
    ----------
    datadir : str
        The path to the directory containing the ensemble-level
        products.
    userdir : str
        The path to the directory containing the user's decisions.
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

    # Setup for user information.
    #
    for key in ['usernotes', 'useraction', 'lastmodified']:
        if not setup_user_setting(jcts, key):
            return None

    # Now add in any user information
    #
    infile = os.path.join(userdir, ensemble,
                          'hull.{}.{:03d}.v{}.json'.format(ensemble,
                                                           mid,
                                                           revision))

    # This is an optional file, so avoid warning messages in the log
    # if we can help it.
    if not os.path.exists(infile):
        return jcts

    ucts = read_json(infile)
    if ucts is None:
        return jcts

    for key in ['usernotes', 'useraction', 'lastmodified']:
        if key in ucts:
            jcts[key]['user'] = ucts[key]

    return jcts


def read_ensemble_json(datadir, userdir, ensemble):
    """Return JSON data from the ensemble.

    The returned structure contains all the versions for this
    ensemble. The "current" version is taken to be the highest
    version number.

    Parameters
    ----------
    datadir : str
        The path to the directory containing the ensemble-level
        products. These are the "proposed" products.
    userdir : str
        The path to the directory containing the user's decisions.
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
                  ['usernotes'] =>
                  ['useraction'] =>
                  ['lastmodified'] =>

           ['usernotes']
           ['lastmodified']

    The usernotes, useraction, lastmodified, and status fields are
    dicts with two keys:
      proposed
      user
    which allows the user to "over ride" the proposed value.

hmmm, the JSON data has {"status": "todo", "stackmap": {"acisfJ1705367m403832_001": 3, "acisfJ1704041m414416_001": 1, "acisfJ1702545m412821_001": 0, "acisfJ1705559m410515_001": 4, "acisfJ1704448m410953_001": 2}, "nmasters": 1, "name": "ens0000900_001", "lastmodified": "", "ncpts": 2, "usernotes": "", "nstacks": 2, "revision": "001"}


    """

    # Process the proposed settings first
    #
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
            hull = read_ensemble_hull_json(datadir, userdir,
                                           ensemble, mid, v)
            if hull is None:
                # if there is a problem reading in a single hull, then
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

        # Override those keys that need proposed/user versions.
        #
        for key in ['usernotes', 'lastmodified', 'status']:
            if not setup_user_setting(jcts, key):
                return None

    revs = list(store['versions'].keys())
    if len(revs) == 0:
        errlog("no JSON data read from ensemble: " +
               "{} {}".format(datadir, ensemble))
        return None

    revs = sorted(revs, key=int, reverse=True)
    store['latest_version'] = revs[0]

    # Now check for user overrides: at present only at the
    # field level.
    #
    inpat = os.path.join(userdir, ensemble, pat)
    matches = glob.glob(inpat)
    if len(matches) == 0:
        return store

    # Assume very-limited metadata here.
    #
    for match in matches:
        jcts = read_json(match)

        revision = jcts['revision']

        try:
            base = store['versions'][revision]
        except KeyError:
            warn("{} has invalid version {}".format(match,
                                                    revision))
            continue

        for k in ['lastmodified', 'usernotes', 'status']:
            try:
                base[k]['user'] = jcts[k]
            except KeyError:
                # for now do not require the user to have all fields
                # set.
                pass

    return store


def read_ensemble_status(datadir, userdir, ensemble, revision):
    """What is the ensemble status?

    What is the current status of the ensemble? This is based on
    read_ensemble_json.


    Parameters
    ----------
    datadir : str
        The path to the directory containing the ensemble-level
        products. These are the "proposed" products.
    userdir : str
        The path to the directory containing the user's decisions.
    ensemble : str
        The ensemble.
    revision : str
        The revision value (in 3-digit form).

    Returns
    -------
    status : {'todo', 'review', 'done'}

    """

    # Try the user and then proposed settings.
    #
    pat = "field.{}.v{}.json".format(ensemble, revision)

    # pick an out-of-bounds value; normally I'd use None but
    # this is a possible value, so use a numeric value.
    #
    not_present = -1

    def lookin(indir):

        infile = os.path.join(indir, ensemble, pat)
        if not os.path.exists(infile):
            return not_present

        try:
            jcts = read_json(infile)
        except IOError:
            return not_present

        try:
            return jcts['status']
        except KeyError:
            return not_present

    ans = lookin(userdir)
    if ans is not None and ans != not_present:
        return ans

    ans = lookin(datadir)
    if ans is not None and ans != not_present:
        return ans

    if ans is None:
        warn("status=None for ensemble " +
             "{} version {}".format(ensemble, revision))
        return "unknown"

    errlog("no status for ensemble " +
           "{} version {}".format(ensemble, revision))
    return "unknown"


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
        errlog("no ensemble dirs found in {}".format(datadir))
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
        jcts = read_ensemble_json(ensembledir + "/../",
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
           'errors': [],
           'usernotes': '',
           'lastmodified': ''}
    for ens in state:
        if 'error' in ens:
            out['errors'].append(ens['ensemble'])
            continue

        if ens['status']['user'] == 'completed':
            key = 'completed'
        elif ens['status']['user'] == 'review':
            key = 'reviews'
        else:
            key = 'todos'

        out[key].append(ens)

    # Has the user saved any notes?
    infile = os.path.join(userdir, 'summary.json')
    if os.path.isfile(infile):
        jcts = read_json(infile)
        if jcts is not None:
            out['usernotes'] = jcts['usernotes']
            out['lastmodified'] = jcts['lastmodified']

    # How about the datatable settings?
    infile = os.path.join(userdir, 'datatable.json')
    if os.path.isfile(infile):
        jcts = read_json(infile)
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


def save_datatable(userdir, data):
    """Save the data table details.

    This could be included in the summary page, but for now
    separate it out.

    Parameters
    ----------
    userdir : str
        The location for the user-stored data.
    data : dict
        The JSON dictionary containing the elements to write out.
    """

    outname = 'datatable.json'
    outfile = os.path.join(userdir, outname)

    with open(outfile, 'w') as fh:
        fh.write(json.dumps(data))


def save_summary(userdir, data):
    """Save the summary details.

    Parameters
    ----------
    userdir : str
        The location for the user-stored data.
    data : dict
        The JSON dictionary containing the elements to write out.
    """

    outname = 'summary.json'
    outfile = os.path.join(userdir, outname)

    store = {"lastmodified": time.asctime(),
             "usernotes": data['usernotes']}
    with open(outfile, 'w') as fh:
        fh.write(json.dumps(store))


def save_ensemble(userdir, data):
    """Save the ensemble-level details.

    Parameters
    ----------
    userdir : str
        The location for the user-stored data.
    data : dict
        The JSON dictionary containing the elements to write out.
    """

    ensemble = data['name']
    version = data['revision']  # this is in string form, 0 padded

    outdir = os.path.join(userdir, ensemble)
    if os.path.exists(outdir):
        if not os.path.isdir(outdir):
            raise IOError("Exists but not a directory! {}".format(outdir))
    else:
        os.mkdir(outdir)

    outname = 'field.{}.v{}.json'.format(ensemble, version)
    outfile = os.path.join(outdir, outname)

    store = {"name": ensemble,
             "lastmodified": time.asctime(),
             "usernotes": data['usernotes'],
             "status": data['status'],
             "revision": version}
    with open(outfile, 'w') as fh:
        fh.write(json.dumps(store))


def save_master(userdir, data):
    """Save the master-level details.

    Parameters
    ----------
    userdir : str
        The location for the user-stored data.
    data : dict
        The JSON dictionary containing the elements to write out.
    """

    ensemble = data['ensemble']
    version = data['revision']  # this is in string form, 0 padded
    masterid = data['masterid']

    outdir = os.path.join(userdir, ensemble)
    if os.path.exists(outdir):
        if not os.path.isdir(outdir):
            raise IOError("Exists but not a directory! {}".format(outdir))
    else:
        os.mkdir(outdir)

    outname = 'hull.{}.{:03d}.v{}.json'.format(ensemble,
                                               masterid,
                                               version)
    outfile = os.path.join(outdir, outname)

    store = {"ensemble": ensemble,
             "masterid": masterid,
             "lastmodified": time.asctime(),
             "usernotes": data['usernotes'],
             "useraction": data['useraction'],
             "revision": version}
    with open(outfile, 'w') as fh:
        fh.write(json.dumps(store))


def save_master_poly(userdir, data):
    """Save the masterhull polygon(s).

    Parameters
    ----------
    userdir : str
        The location for the user-stored data.
    data : dict
        The JSON dictionary containing the elements to write out.
    """

    ensemble = data['ensemble']
    version = data['revision']  # this is in string form, 0 padded
    masterid = data['masterid']

    outdir = os.path.join(userdir, ensemble)
    if os.path.exists(outdir):
        if not os.path.isdir(outdir):
            raise IOError("Exists but not a directory! {}".format(outdir))
    else:
        os.mkdir(outdir)

    outname = 'poly.{}.{:03d}.v{}.json'.format(ensemble,
                                               masterid,
                                               version)
    outfile = os.path.join(outdir, outname)

    store = {"ensemble": ensemble,
             "masterid": masterid,
             "lastmodified": time.asctime(),
             "polygons": data['polygons'],
             "revision": version}
    with open(outfile, 'w') as fh:
        fh.write(json.dumps(store))


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

    state = read_ensemble_json(datadir, userdir, ensemble)
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
        out += "- see Doug!</p></body></html>"
        return 404, out

    hullfile = os.path.join(dname,
                            'master_hulls.' +
                            '{}.v{}.fits'.format(ensemble, revstr))
    if not os.path.isfile(hullfile):
        errlog("missing hullfile {}".format(hullfile))
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
    # NOTE: there is code that does something similar to this
    #       above
    #
    try:
        ds = pycrates.CrateDataset(hullfile, mode='r')
    except IOError as exc:
        errlog("unable to read hullfile {} - {}".format(hullfile, exc))
        out = "<!DOCTYPE html><html><head><title>INTERNAL ERROR</title>"
        out += "</head><body><p>Error reading "
        out += "hullfile={}\nreason=\n{}".format(hullfile, exc)
        out += "- see Doug!</p></body></html>"
        return 404, out

    try:
        cr = ds.get_crate('HULLMATCH')
    except IndexError as exc:
        errlog("unable to read HULLMATCH from " +
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
        errlog("missing STKIDNUM keyword in hullfile {} - {}".format(hullfile))
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
            errlog("missing {} keyword in hullfile {} - {}".format(key, hullfile))
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
            errlog("unexpected stackid={}".format(stackid))
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
        errlog("unable to read HULLLIST from " +
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

            wcs_orig = []
            for nvertex, eqpos in zip(qcr.NVERTEX.values,
                                      qcr.EQPOS.values):
                wcs_orig.append(eqpos_to_dict(eqpos, nvertex))

            qcr = None

        store['wcs_orig'] = wcs_orig

        # Has the user got their own version?
        polyfile = os.path.join(userdir, ensemble,
                                'poly.{}.{:03d}.v{}.json'.format(ensemble,
                                                                 midval,
                                                                 revstr))
        if os.path.exists(polyfile):
            jcts = read_json(polyfile)
            if jcts is None:
                errlog("polyfile is unreadable: {}".format(polyfile))
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
            log("missing stack region file {}".format(regfile))
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
            log("key is missing in mancode_by_component: {}".format(key))
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
            errlog("Unable to open xmdat file: {}".format(infile))
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

        components.append({'name': name,
                           'eband': ebands_by_component[key],
                           'likelihood': likelihood_by_component[key],
                           'adjusted': mancode_by_component[key] > 0,
                           'svdqa': svdqa_by_component[key],
                           'include_in_centroid': centroid_by_component[key],
                           'mrg3rev': mrgrev_by_component[key]})

    # Do we know the actual size here?
    #
    if len(components) == 0:
        log("filtered out all components!")
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
        userdir = context['userdir']
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
                xmdatdir = context['xmdatdir']
                status, cts = create_master_hull_page(env,
                                                      datadir,
                                                      userdir,
                                                      datadir,
                                                      xmdatdir,
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

        upath = urlparse.urlparse(self.path)
        path = upath.path

        # Strip the leading /
        if path != "":
            path = path[1:]

        try:
            savefunc = {'save/summary': save_summary,
                        'save/ensemble': save_ensemble,
                        'save/master': save_master,
                        'save/masterpoly': save_master_poly,
                        'save/datatable': save_datatable
                        }[path]
        except KeyError:
            errlog("Unexpected POST path={}".format(path))
            self.send_error(404)
            return

        # for now require JSON
        ctype = self.headers['Content-Type'].split(';')
        if not ctype[0] == 'application/json':
            errlog("Unexpected content-type: {}".format(ctype[0]))
            self.send_error(404)
            return

        clen = int(self.headers['Content-Length'])
        data = self.rfile.read(clen)
        try:
            jcts = json.loads(data)
        except Exception as exc:
            errlog("Invalid JSON: {}".format(exc))
            self.send_error(404)
            return

        context = self.server.context
        userdir = context['userdir']

        try:
            savefunc(userdir, jcts)
        except Exception as exc:
            errlog("Unable to save data: {}".format(exc))
            self.send_error(404)
            return

        self._set_headers(200)
        self.wfile.write("")

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
            cts = read_ensemble_json(datadir, userdir, ensemble)
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
            errlog("Invalid master id: {}".format(toks[2]))
            self.send_error(404)

        ensemble_status = read_ensemble_status(datadir, userdir,
                                               ensemble, revision)
        cts = read_ensemble_hull_json(datadir, userdir,
                                      ensemble, masterid, revision)
        if cts is None:
            self.send_error(404)
        else:
            cts['ensemble_status'] = ensemble_status
            self.send_as_json(cts)

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
            errlog("failed glob={}".format(pat))
            self.send_error(404)
            return

        # Should pick the highest version; for now pick the first
        #
        if len(matches) > 1:
            log("Multiple matches for {}".format(pat))

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
          xmdatdir='/data/L3/chs_master_match/input/xmdat3',
          evt3dir='/data/L3/chs_master_match/input/stkevt3'):

    for dirname in [userdir, webdir, datadir, evt3dir, xmdatdir,
                    templatedir]:
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
    log("Web assets: {}".format(webdir))
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

    datadir = os.path.normpath(os.path.abspath(sys.argv[1]))

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
    # within the data directory, to avois over-writing files.
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
          port=port)
