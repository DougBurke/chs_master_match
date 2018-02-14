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

try:
    import pycrates
except ImportError:
    sys.stderr.write("Needs crates!\n")
    sys.exit(1)


def log(msg):
    """Log a message.

    Parameters
    ----------
    msg: str
        The message
    """

    print("LOG: {}".format(msg))


def errlog(msg):
    """Log an error message.

    Parameters
    ----------
    msg: str
        The message
    """

    log("ERROR: {}".format(msg))


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
        print("LOG: error reading JSON from {}\n{}".format(infile, exc))
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
        print("LOG: no field.json files found for ensemble " +
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
            print("LOG: missing revision keyword in {}".format(match))
            continue

        # this should not happen, so do not worry too much about the
        # error handler
        if v in store['versions']:
            print("LOG: multiple revision={} in {}".format(v, match))
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
            print("LOG: no masters for revsision {} in {}".format(v, match))
            continue

        if 'masters' in jcts:
            print("LOG: overwriting masters setting in {}".format(match))

        jcts['masters'] = hulls
        store['versions'][v] = jcts

    revs = list(store['versions'].keys())
    if len(revs) == 0:
        print("LOG: no JSON data read from ensemble: {} {}".format(datadir, ensemble))
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
        log("LOG: no ensemble dirs found in {}".format(datadir))
        return None

    # HACK: use first 20 ensembles to save time
    #
    log("DEBUG: restricting to first 20 ensembles")
    ensembledirs = ensembledirs[:20]

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
                        print("LOG: expected coordsys, found " +
                              "[{}] in {}".format(l, infile))
                        return None

                    csys = l
                    continue

                if not l.startswith('polygon('):
                    print("LOG: expected polygon(...), found " +
                          "[{}] in {}".format(l, infile))
                    return None

                return "{}; {}".format(csys, l)

    except IOError as exc:
        print("LOG: can not read region file:" +
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
        print("ERROR: no row matching {}".format(fname))
        return None

    elif nrows > 1:
        print("ERROR: multiple rows matching {}".format(fname))
        return None

    # NOTE: QA cases have NVERTEX=0
    nv = cr.NVERTEX.values[0]
    if nv < 3:
        print("ERROR: nvertex={} for {}".format(nv, fname))
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
        print("LOG: ensemble json missing key {}".format(exc))
        return None

    mid = "{:03d}".format(masterid)
    hull = None
    for m in masters:
        if m['masterid'] == mid:
            hull = m
            break

    if hull is None:
        print("LOG: missing hull {}".format(mid))
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
        print("LOG: ERROR missing rawdir {}".format(dname))
        return None

    hullfile = os.path.join(dname,
                            'master_hulls.' +
                            '{}.v{}.fits'.format(ensemble, version))
    if not os.path.isfile(hullfile):
        print("LOG: ERROR missing hullfile {}".format(hullfile))
        return None

    infile = "{}[master_id={}][cols stackid, component]".format(hullfile,
                                                                masterid)
    cr = pycrates.read_file(infile)
    nrows = cr.get_nrows()
    if nrows == 0:
        print("LOG: ERROR no hull data from {}".format(infile))
        return None

    # Just warn, nbut do not error out, if there's a difference
    #
    if nrows != hull['ncpts']:
        print("LOG: WARNING expected {} but got {} from {}".format(hull['ncpts'],
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


def create_index_page(datadir):
    """Create the top-level page: ensemble status.

    Parameters
    ----------
    datadir : str
        The path to the directory containing the ensemble-level
        products.

    Returns
    -------
    status, page : int, str
        The response status and HTML contents
    """

    out = "<DOCTYPE html><html><head><meta charset='UTF-8'>"
    out += "<title>CHS review</title>"
    out += """<style type='text/css'>
.selection {
  display: grid;
  grid-gap: 0.2em;
  grid-template-columns: repeat(auto-fill, 13em);

}

.label {
  font-size: larger;
}

.selection {
  border: 1px solid black;
  padding: 0.2em 0.5em;
  min-height: 3em;
}

.ensemble {
  border: 1px solid black;
  display: inline-block;
  padding: 0.2em 0.5em;
  text-decoration: none;
}

.ensemble:hover {
  background-color: rgba(0, 0, 0, 0.2);
  cursor: pointer;
}

#summary {
  display: grid;
  grid-template-columns: 7em 1fr;
}
</style><script type='text/javascript'>
function plural(x) { if (x === 1) { return ""; } else { return "s"; } }

function updatePage(json) {
  var ntodos = json.todos.length;
  var nreviews = json.reviews.length;
  var ncompleted = json.completed.length;

  function nens(n) { return n.toString() + " ensemble" + plural(n); }
  document.getElementById("ntodos").innerHTML = nens(ntodos);
  document.getElementById("nreviews").innerHTML = nens(nreviews);
  document.getElementById("ncompleted").innerHTML = nens(ncompleted);

  // TODO: do we have to clear out these divs?
  var parent = document.getElementById("todo");
  for (let i = 0; i < ntodos; i++) {
    var ens = json.todos[i];

    // TODO: need to work out the number of remaining hulls
    //       or just do away with this concept (but would be nice
    //       to know which ones you've worked on)
    //
    var el = document.createElement("a");
    el.className = "ensemble";
    el.href = ens['name'];
    el.innerHTML = ens['name'] + "<br>" +
        ens['nmasters'].toString() +
        " remaining of " +
        ens['nmasters'].toString() + ", " +
        ens['nstacks'].toString() + " stack" + plural(ens['nstacks']) +
        "<br>N/A";

    parent.appendChild(el);
  }

  parent = document.getElementById("review");
  for (let i = 0; i < nreviews; i++) {
    var ens = json.reviews[i];

    // TODO: need to work out the number of remaining hulls
    var el = document.createElement("a");
    el.className = "ensemble";
    el.href = ens['name'];
    el.innerHTML = ens['name'] + "<br>" +
        ens['nmasters'].toString() + " hulls, " +
        ens['nstacks'].toString() + " stack" + plural(ens['nstacks']);

    parent.appendChild(el);
  }

  parent = document.getElementById("completed");
  for (let i = 0; i < ncompleted; i++) {
    var ens = json.completed[i];

    // TODO: need to work out the number of remaining hulls
    var el = document.createElement("a");
    el.className = "ensemble";
    el.href = ens['name'];
    el.innerHTML = ens['name'] + "<br>" +
        ens['nmasters'].toString() + " hulls, " +
        ens['nstacks'].toString() + " stack" + plural(ens['nstacks']);

    parent.appendChild(el);
  }
}

function initialize() {
  var httpRequest = new XMLHttpRequest();
  if (!httpRequest) {
      alert("Unable to create a XMLHttpRequest!");
      return;
  }
  httpRequest.addEventListener("load", function() {
    updatePage(httpRequest.response);
  });
  httpRequest.addEventListener("error", function() {
      alert("Unable to load data!");
  });
  httpRequest.addEventListener("abort", function() {
      alert("Unable to load data!");
  });

  // Do I need to add a cache-busting identifier?
  httpRequest.open('GET', '/api/summary?' +
                   (new Date()).getTime());
  httpRequest.responseType = 'json';
  httpRequest.send();
}
</script></head><body onload='initialize();'>"""

    out += "<h1>Ensemble selection page</h1>"
    out += "<div id='summary'>"
    # TODO: update this value
    out += "<span>Last save:</span><span>never</span>"
    out += "<span>To do:</span>"
    out += "<span id='ntodos'></span>"
    out += "<span>To review:</span>"
    out += "<span id='nreviews'></span>"
    out += "<span>Completed:</span>"
    out += "<span id='ncompleted'></span>"
    out += "<span>Directory:</span><span>{}</span>".format(datadir)
    out += "</div>"

    # notes
    out += "<p class='label'>Notes</p><textarea cols=80 rows=10 "
    out += "id='usercontent'></textarea><br><button "
    out += "id='saveusercontent' class='button'>Save notes</button>"

    # To do
    #
    # Could allow sorting by several different metrics; use
    # data elements to encode the values for the sort?
    #
    out += "<p class='label'>To do</p><div id='todo' class='selection'>"
    out += "</div>"

    # Review
    out += "<p class='label'>Review</p><div id='review' class='selection'>"
    out += "</div>"

    # Completed
    out += "<p class='label'>Completed</p><div id='completed' class='selection'>"
    out += "</div>"

    out += "</body></html>"
    return 200, out


def create_ensemble_index_page(datadir, ensemble):
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
        out = "<DOCTYPE html><html><head><title>ERROR</title>"
        out += "</head><body><p>There was an error when processing "
        out += "the directory {} for ".format(datadir)
        out += "ensemble {}".format(ensemble)
        out += "</body></html>"
        return 404, out

    out = "<DOCTYPE html><html><head><meta charset='UTF-8'>"
    out += "<title>{}</title>".format(ensemble)
    out += """<style type='text/css'>
#infobar {
  float: left;
  width: 20em;
}

#infobar div {
  border: 1px solid black;
  padding: 0.2em 0.5em;
  margin-bottom: 1em;
}

#overviewbar {
  float: left;
}

.info {
  display: grid;
  grid-gap: 0.2em;
  grid-template-columns: 7em 7em;

}

.label {
  grid-column: 1 / 3;
  /* text-align: center; */
}

#username {
  border: 1px solid black;
  float: right;
  font-size: x-large;
  padding: 0.2em 0.5em;
}

th {
  background: rgba(102, 204, 102, 0.4);
}

tr:nth-child(even) {
  background: rgba(102, 102, 204, 0.4);
}

/* Ensure the decision field fills all the space */
td:nth-child(4) {
  width: 100%;
}

td.undecided {
  background-color: orange;
}

#final {
  text-align: center;
}

#finish {
}

#home {
  text-align: center;
}

#home a {
  display: inline-block;
  padding: 0.2em 0.5em;
}

#home a:hover {
}
</style><script type='text/javascript'>
var state;

function removeChildren(parent) {
  while (parent.firstChild) {
    parent.removeChild(parent.firstChild);
  }
}

function updatePage(json) {
  state = json;

  // Set the version options for the user to select.
  //
  var versions = [];
  for (let revision in state.versions) {
    versions.push(revision);
  }

  // wasteful comparison search but should not be an issue here
  versions.sort((a, b) => {
    const aa = Math.floor(a);
    const bb = Math.floor(b);
    if (aa < bb) { return -1; }
    if (aa > bb) { return 1; }
    return 0;
    }).reverse();

  var parent = document.getElementById('version');
  removeChildren(parent);
  for (let i = 0; i < versions.length; i++) {
    const el = document.createElement('option');
    el.value = versions[i];
    if (i === 0) { el.selected = true; }
    el.innerHTML = versions[i];

    parent.appendChild(el);
  }

  if (versions.length === 1) {
    parent.disabled = true;
  } else {
    parent.addEventListener("change", (e) => { setVersion(e.target.value); });
  }

  // default to the latest version
  setVersion(state.latest_version);
}

// Switch the display to the given version; minimal error checking
//
// If not the latest then page is "read-only".
//
function setVersion(revision) {

  const isLatest = state.latest_version === revision;

  let info = state.versions[revision];
  if (typeof info === 'undefined') {
    alert("Internal error: unknown revision: [" + revision + "]");
    return;
  }

  document.getElementById("nstacks").innerHTML = info.nstacks.toString();
  document.getElementById("nmasters").innerHTML = info.nmasters.toString();

  /* handle the master list */
  let parent = document.getElementById('summarydata');
  let all_hulls_have_useraction = true;
  removeChildren(parent);
  for (let i = 0; i < info.nmasters; i++) {
    let m = info.masters[i];
    let el = document.createElement("tr");

    var td;
    td = document.createElement("td");
    td.innerHTML = m.masterid;
    el.appendChild(td);

    td = document.createElement("td");
    td.innerHTML = m.ncpts.toString();
    el.appendChild(td);

    td = document.createElement("td");
    var a = document.createElement("a");
    a.className = 'hullreview';
    a.href = "/" + info.name + "/" + revision + "/" + m.masterid;
    if (isLatest) { a.innerHTML = "Review"; } else { a.innerHTML = "View"; }
    td.appendChild(a);
    el.appendChild(td);

    td = document.createElement("td");
    if (m.useraction === '') {
      td.className = "undecided";
      td.innerHTML = "NO DECISION";
      all_hulls_have_user_action = false;
    } else {
      td.innerHTML = m.useraction;
    }
    el.appendChild(td);

    parent.appendChild(el);
  }

  // Update the image
  parent = document.getElementById('overviewimage');
  parent.src = "/img/" + info.name + "/field." + info.name +
               ".v" + revision + ".png";

  // Is this an actionable page or not?
  // a) are we the latest version
  // b) do all the hulls have a user action?
  //
  const finish = document.getElementById('final');
  if (isLatest) {
    finish.style.display = 'block';
    document.getElementById('finish')
            .disabled = !all_hulls_have_user_action;;
  } else {
    finish.style.disply = 'none';
  }
}

function initialize() {
  let httpRequest = new XMLHttpRequest();
  if (!httpRequest) {
      alert("Unable to create a XMLHttpRequest!");
      return;
  }
  httpRequest.addEventListener("load", function() {
    updatePage(httpRequest.response);
  });
  httpRequest.addEventListener("error", function() {
      alert("Unable to load data!");
  });
  httpRequest.addEventListener("abort", function() {
      alert("Unable to load data!");
  });

  // Do I need to add a cache-busting identifier?
  httpRequest.open('GET', '/api/""" + ensemble + """?' +
                   (new Date()).getTime());
  httpRequest.responseType = 'json';
  httpRequest.send();
}
</script></head><body onload='initialize();'>"""

    out += "<div id='infobar'>"

    out += "<div id='home'><a href='/'>Ensemble list</a></div>"

    out += "<div class='info'>"
    out += "<span class='label'>{}</span>".format(ensemble)
    out += "<span>Version:</span>"

    # If only one version then could not display a menu, but
    # that complicates the code a bit so leave for now
    # (may need to change given user feedback)
    #
    out += "<select id='version' class='button'>"
    out += "</select>"

    out += "<span>Stacks:</span>"
    out += "<span id='nstacks'></span>"

    out += "<span>Master Hulls:</span>"
    out += "<span id='nmasters'></span>"

    out += "</div>"

    out += "<div id='summary'><table><thead><tr>"
    for n in ['Master', 'Ncpt', 'Action', 'Decision']:
        out += "<th>{}</th>".format(n)

    out += "</tr></thead><tbody id='summarydata'>"
    out += "</tbody></table></div>"

    # notes
    out += "<div id='notes'><p class='label'>Notes</p>"
    out += "<textarea cols=32 rows=10 id='usercontent'></textarea>"
    out += "<br><button id='saveusercontent' class='button'>"
    out += "Save notes</button></div>"

    # final
    out += "<div id='final'><button id='finish' class='button' "
    out += "type='button'>Finish review</button></div>"

    # overview image
    out += "</div><div id='overviewbar'>"
    out += "<div id='username'>{}</div>".format(os.getlogin())
    out += "<img id='overviewimage'>"
    out += "</div>"

    out += "</body></html>"
    return 200, out


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


def create_master_hull_page(datadir, rawdir,
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
        out = "<DOCTYPE html><html><head><title>ERROR</title>"
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
        out = "<DOCTYPE html><html><head><title>ERROR</title>"
        out += "</head><body><p>Invalid master hull.</p></body></html>"
        return 404, out

    elif len(hulls) > 1:
        print("LOG: hulls = {}".format(hulls))
        out = "<DOCTYPE html><html><head><title>INTERNAL ERROR</title>"
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
        print("LOG: ERROR missing rawdir {}".format(dname))
        out = "<DOCTYPE html><html><head><title>INTERNAL ERROR</title>"
        out += "</head><body><p>Missing dir={}".format(dname)
        out += + "- see Doug!</p></body></html>"
        return 404, out

    hullfile = os.path.join(dname,
                            'master_hulls.' +
                            '{}.v{}.fits'.format(ensemble, revstr))
    if not os.path.isfile(hullfile):
        print("LOG: ERROR missing hullfile {}".format(hullfile))
        out = "<DOCTYPE html><html><head><title>INTERNAL ERROR</title>"
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
        print("LOG: unable to read hullfile {} - {}".format(hullfile, exc))
        out = "<DOCTYPE html><html><head><title>INTERNAL ERROR</title>"
        out += "</head><body><p>Error reading "
        out += "hullfile={}\nreason=\n{}".format(hullfile, exc)
        out += + "- see Doug!</p></body></html>"
        return 404, out

    try:
        cr = ds.get_crate('SRCMATCH')
    except IndexError as exc:
        print("LOG: unable to read SRCMATCH from hullfile {} - {}".format(hullfile, exc))
        out = "<DOCTYPE html><html><head><title>INTERNAL ERROR</title>"
        out += "</head><body><p>Error reading SRCMATCH block of "
        out += "hullfile={}\nreason=\n{}".format(hullfile, exc)
        out += + "- see Doug!</p></body></html>"
        return 404, out

    # We have aleady got this information, but recreate it
    #
    nstks = cr.get_key_value('STKIDNUM')
    if nstks is None:
        print("LOG: missing STKIDNUM keyword in hullfile {} - {}".format(hullfile))
        out = "<DOCTYPE html><html><head><title>INTERNAL ERROR</title>"
        out += "</head><body><p>No STKIDNUM keyword in "
        out += "hullfile={}\n".format(hullfile)
        out += + "- see Doug!</p></body></html>"
        return 404, out

    stack_map = {}
    for i in range(nstks):
        key = "STKID{:03d}".format(i)
        val = cr.get_key_value(key)
        if val is None:
            print("LOG: missing {} keyword in hullfile {} - {}".format(key, hullfile))
            out = "<DOCTYPE html><html><head><title>INTERNAL ERROR</title>"
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
        print("LOG: unable to read SRCLIST from hullfile {} - {}".format(hullfile, exc))
        out = "<DOCTYPE html><html><head><title>INTERNAL ERROR</title>"
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
                print("LOG: ERROR missing qafile {}".format(qafile))
                out = "<DOCTYPE html><html><head><title>INTERNAL ERROR</title>"
                out += "</head><body><p>Missing qafile={}".format(qafile)
                out += "- see Doug!</p></body></html>"
                return 404, out

            try:
                qcr = pycrates.read_file(qafile + "[cols nvertex, eqpos]")
            except IOError as exc:
                print("LOG: ERROR unable to read qafile {}".format(qafile))
                out = "<DOCTYPE html><html><head><title>INTERNAL ERROR</title>"
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
            print("LOG missing stack region file {}".format(regfile))
            out = "<DOCTYPE html><html><head><title>INTERNAL ERROR</title>"
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
                # THIS SHOULD BE CONSIDERED AN ERROR
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
            pass

        store.append(shull)

    out = "<DOCTYPE html><html><head><meta charset='UTF-8'>"
    out += "<title>{}: v{} hull {}</title>".format(ensemble,
                                                   revstr,
                                                   mid)

    # JS9 support
    out += '<link type="text/css" rel="stylesheet" href="/js/js9/js9support.css">'
    out += '<link type="text/css" rel="stylesheet" href="/js/js9/js9.css">'
    out += '<link rel="apple-touch-icon" href="/js/js9/images/js9-apple-touch-icon.png">'
    out += '<script type="text/javascript" src="/js/js9prefs-chs.js"></script>'
    out += '<script type="text/javascript" src="/js/js9/js9support.min.js"></script>'
    out += '<script type="text/javascript" src="/js/js9/js9.min.js"></script>'
    out += '<script type="text/javascript" src="/js/js9/js9plugins.js"></script>'

    # Extra JS code
    #
    # This adds window.convexHull
    #
    out += '<script type="text/javascript" src="/js/convexhull.js"></script>'

    # WebSAMP
    #
    out += '<script type="text/javascript" src="/js/samp.js"></script>'

    # Some of the javascript code - handling multiple pages - could
    # only be included when needed (since it is known at this point
    # whether it is needed), but it's not worth the effort to
    # support this.
    #
    out += """<style type='text/css'>
#infobar {
  float: left;
  width: 20em;
}

#infobar div {
  border: 1px solid black;
  padding: 0.2em 0.5em;
  margin-bottom: 1em;
}

#overviewbar {
  float: left;
}

.info {
  display: grid;
  grid-gap: 0.2em;
  grid-template-columns: 10em 7em;

}

.label {
  grid-column: 1 / 3;
  /* text-align: center; */
}

#username {
  border: 1px solid black;
  float: right;
  font-size: x-large;
  padding: 0.2em 0.5em;
}

th {
  background: rgba(102, 204, 102, 0.4);
}

tr:nth-child(even) {
  background: rgba(102, 102, 204, 0.4);
}

/* Ensure the decision field fills all the space */
td:nth-child(4) {
  width: 100%;
}

td.undecided {
  background-color: orange;
}

#action {
  text-align: center;
}

#summary {
  text-align: center;
}

select + table {
  margin-top: 1em;
}

#usercontent {
  min-height: 3em;
}

#home {
  text-align: center;
}

#home a {
  display: inline-block;
  padding: 0.2em 0.5em;
}

#home a:hover {
}

.stackid {
  font-size: larger;
  min-height: 1.5em;
}

.rebin , .blur {
  background-color: rgba(217, 83, 79, 0.6);  /* red-ish */
  border-radius: 0.3em;
  float: left;
  margin-right: 0.5em;
  padding-bottom: 0.2em;
  padding-left: 0.2em;
  padding-top: 0.2em;
}

.reload , .zoom {
  float: right;
}

/* The "label" for the button bar */
.rebin span , .blur span {
  border-right-color: rgba(20, 20, 20, 0.4);
  border-right-style: solid;
  border-right-width: 3px;
  font-family: monospace;
  font-size: smaller;
  padding-bottom: 0.2em;
  padding-right: 0.2em;
  padding-top: 0.2em;
}

/* restyle the radio buttons to try and save space */
div.useropts input[type='radio'] {
  display: none;
}

div.useropts label {
  border: 1px solid rgba(212, 63, 58, 0.4);
  color: white;
  /* cursor: pointer; */
  /* display: inline-block; */
  font-family: monospace;
  font-size: smaller;
  padding-bottom: 0.2em;
  padding-left: 0.4em;
  padding-right: 0.4em;
  padding-top: 0.2em;
}

div.useropts label:hover {
  background-color: rgba(210, 35, 45, 0.6);
  border-color: rgb(172, 41, 37);
}

div.useropts input[type='radio']:checked + label {
  background-color: rgba(217, 83, 79, 0.8);  /* red-ish */
  border: 1px solid rgb(212, 63, 58);
  color: black;
}

.hulllink {
  border: 1px solid rgb(212, 63, 58);
  border-radius: 0.5em;
  color: white;
  display: inline-block;
  font-family: monospace;
  font-size: smaller;
  margin-top: 0.4em;
  padding: 0.2em 0.5em;
  width: 2em;
}

span.hulllink {
  color: black;
}

a.hulllink {
  background-color: rgba(217, 83, 79, 0.8);  /* red-ish */
  cursor: pointer;
  text-decoration: none;
}

a.hulllink:hover {
  background-color: rgba(210, 35, 45, 0.8);
  border-color: rgb(172, 41, 37);
}

.hulllink + .hulllink {
  margin-left: 0.2em;
}
</style><script type='text/javascript'>
var ensemble = '"""

    out += ensemble
    out += """';
var revision = '"""
    out += revstr
    out += """';
var npages = """
    out += str(hull['npages'])
    out += """;
var state;

var pageNum;
var imageScale = 'log10';

// Change to the given page and scaling.
//
function changePage() {
  let pageStr = pageNum.toString();
  if (pageNum < 10) { pageStr = "00" + pageStr; }
  else if (pageNum < 100) { pageStr = "0" + pageStr; }

  const img = '/img/' + ensemble + '/hull.' + ensemble + '."""

    out += mid + ".p' + pageStr + '.v" + revstr
    out += """.' + imageScale + '.png';

  document.getElementById('pageview').src = img;
}

function setPage(newPage) {
  if ((newPage < 1) || (newPage > npages)) {
    alert("Internal error: newPage = " + newPage.toString());
    return;
  }

  pageNum = newPage;
  changePage();
}

function setScaling(newScale) {
  imageScale = newScale;
  changePage();
}

// Map from stack value to the detected energy band.
// This is tricky because for multiple components there is
// no guarantee that the band is the same for the different
// components. So the choice is either to go for a component
// view, or to either pick one band or combine them.
// Thew assumption is that this has already been decided by
// the time this JS code has been written out.
//
var enbands = """

    out += json.dumps(stack_band)
    out += """;

// Store the region information for this *ensemble*; that is,
// provide coordinates (WCS) for all master hulls and stacks.
//
// This allows the UI to then decide whether to show all this
// information or not.
//
var regionstore = {"""

    out += "'masterhulls': " + json.dumps(hull_store) + ","
    out += "'stackhulls': " + json.dumps(stack_polys)
    out += """};

// pan to this position; opts is the argument to pass
// to JS9 commands to determine the window to use.
// wcs contains ra and dec fields.
//
function goToRaDec(wcs, opts) {
  const pix = JS9.WCSToPix(wcs.ra, wcs.dec, opts);
  JS9.SetPan(pix.x, pix.y, opts);
}

// Finalize the widgets in the window and add the regions.
//
function finalizeJS9Display(ensemble, masterid, stack, winid) {

  return function(img) {
    setupJS9(img, ensemble, masterid, stack, winid);
    addRegionsToJS9(img, ensemble, masterid, stack, winid);

/***
    const out = addRegionToJS9(img, ensemble, masterid, stack, winid,
                             data);
    updateJS9StackCounter(winid, out.nhulls);
    const opts = {display: img};
    if (typeof(out.center) !== "undefined") {
      goToRaDec(out.center, opts);
    }
***/
  };
}

function _old_finalizeJS9Display(ensemble, masterid, stack, winid) {
  return function (img) {
    $.ajax({url: '/regions/ensemble/' + ensemble +
                 '/' + masterid.toString(),
            dataType: 'json'})
     .done((data, textStatus) => {
       setupJS9Widgets(img, ensemble, masterid, stack, winid);
       const out = addRegionToJS9(img, ensemble, masterid, stack, winid,
                                  data);
       updateJS9StackCounter(winid, out.nhulls);
       const opts = {display: img};
       if (typeof(out.center) !== "undefined") {
          goToRaDec(out.center, opts);
       }
     });
  };
}

//   - stack      for stack-level hulls
//   - original   the input 'master-level' hulls
//   - convex     the convex hull calculated from the current hull
//   - regions    the user-editable hull
//                i.e. this is the default regions layer
//
// Only the latter is user-editable or selectable (at least
// directly). For now we use the "standard" regions layer
// for this
//
const stackLayer = 'stack';
const convexLayer = 'convex';
const originalLayer = 'original';
const masterLayer = 'regions';

// Customize the JS9 display window
//
// winid is the base HTML id of the "light window" containing
// the JS9 display.
//
// Is also sets up the layers.
//
function setupJS9(img, ensemble, masterid, stack, winid) {

  // Hard code logic for the parent element (could find it by
  // hunting up the tree but can not be bothered at the moment)
  //
  const container = document.getElementById('d' + winid);

  // Add handlers for the user buttons:
  //   - blur
  //   - rebin
  //   - reload region
  //   - zoom
  //
  const opts = {display: img};
  let btns = container.getElementsByClassName("sigma");
  for (let i = 0; i < btns.length; i++) {
    btns[i].addEventListener("change", (e) => {
      JS9.GaussBlurData(e.target.value, opts);
    });
  }

  const sigma0 = document.getElementById(winid + "sigma0");

  btns = container.getElementsByClassName("binsize");
  for (let i = 0; i < btns.length; i++) {
    btns[i].addEventListener("change", (e) => {

      // What is the current position - I am not too bothered
      // if I do not calculate the exact center, as there is
      // an effective zoom in/out going on here due to the
      // change in binning - but let's see if this is visually
      // confusing (the issue is if I am really calculating the
      // center of the image here).
      //
      const idata = JS9.GetImageData(false, opts);

      // This does not seem to be working as intended!
      const wcs = JS9.PixToWCS(idata.width / 2, idata.height / 2, opts);

      JS9.DisplaySection({bin: e.target.value}, opts);

      // Could perhaps just jump to the master hull
      // goToRaDec(wcs, opts);  not working well

      // reset the blur button to 0 since rebinning removes
      // the blurring automatically
      sigma0.checked = true;
    });
  }

  // When reloading the regions, delete all the existing ones
  // as a precaution. The region info is reloaded from disk
  // (an alternative would be to just cache the data and re-use
  // it).
  //
  document.getElementById(winid + 'ReloadRegions')
    .addEventListener("click", (e) => {
      JS9.RemoveRegions("all", opts);
      $.ajax({url: '/regions/ensemble/' + ensemble +
                   '/' + masterid.toString(),
              dataType: 'json'})
       .done((data, textStatus) => {
         addRegionToJS9(img, ensemble, masterid, stack, winid, data);
       });
    });

  container.getElementsByClassName("zoomin")[0]
     .addEventListener("click", (e) => { JS9.SetZoom("in", opts); });
  container.getElementsByClassName("zoomout")[0]
     .addEventListener("click", (e) => { JS9.SetZoom("out", opts); });

  // toggle the panner; it looks like need to give winid to display,
  // not img.
  //
  document.getElementById(winid + 'ShowPanner')
      .addEventListener("click", (e) => { JS9.DisplayPlugin('panner',
                                                            {display: winid}); });

  // Set up the layers. We have to base the options on
  // JS9.Regions.opts to get "sensible" behavior.
  //
  let layerOpts = Object.assign({}, JS9.Regions.opts);
  for (let name of ['movable', 'rotatable', 'resizable', 'evented']) {
    layerOpts[name] = false;
  }

  // Make sure we remove the onchange handler from these extra layers
  layerOpts.onchange = null;

  for (let name of [stackLayer, originalLayer, convexLayer]) {
    JS9.NewShapeLayer(name, layerOpts, opts);
  }
}

// The layer argument was an attempt to allow shapes in different
// layers, to better support interactivity, but it doesn't quite
// work right as is and I have a (hopefully) working solution
// for the interactivity.
//
function add_hull_to_js9(hull, opts, win, layer='regions') {

  // Need to convert to image coordinates
  const pts = [];
  for (let i = 0; i < hull.ra.length; i++) {
    pts.push(JS9.WCSToPix(hull.ra[i], hull.dec[i], win));
  }

  // Need to set the label
  const shape = {shape: 'polygon', pts: pts};
  if (typeof(hull.label) !== "undefined") {
    shape.text = hull.label;
  }
  JS9.AddShapes(layer, shape, opts, win);
  // JS9.AddRegions(shape, opts, win);

}


// Add the stack-level and master hull(s) to the JS9 window.
//
// Stack level hulls are drawn first (maybe in a different layer?)
//
function addRegionsToJS9(img, ensemble, masterid, stack, regions) {

  let stackhulls = regionstore.stackhulls[stack];
  if (typeof stackhulls === "undefined") {
    alert("no stack-level hulls found for " + stack);
    return;
  }

  let display = {display: img};

  let linestyle = [1];
  let hullOpts = {color: 'orange',
                  strokeDashArray: linestyle,
                  changeable: false,
                  tags: 'stack'};

  for (let shull of stackhulls) {
    // linestyle: solid for mancode is 0, otherwise
    // dotted.
    //
    // wanted to add strokeWidth to increase the width, but can
    // not get it to work
    //
    if (shull.mancode === 0) {
      delete hullOpts.strokeDashArray;
    } else {
      hullOpts.strokeDashArray = linestyle;
    }
    add_hull_to_js9(shull, hullOpts, display, stackLayer);
  }

  /*
   * TODO: the tag should act as a unique identifier
   *       which will be useful when supporting QA cases
   *       (when can have multiple hulls).
   */
  let masterhull = regionstore.masterhulls[masterid];
  if (typeof masterhull === "undefined") {
    alert("No master hull found!");
    return;
  }

  if (masterhull.wcs.length === 0) {
    console.log("We have no master hull information for " +
                ensemble + " " + masterid.toString() + "!");
    return;
  }

  hullOpts = {movable: false,
              rotatable: false,
              resizable: false,
              tags: 'master'};

  let ras = [];
  let decs = []

  const origOpts = Object.assign({}, hullOpts);
  origOpts.color = 'white';
  origOpts.strokeDashArray = [3, 3];

  for (let hull of masterhull.wcs) {
    add_hull_to_js9(hull, origOpts, display, originalLayer);
    add_hull_to_js9(hull, hullOpts, display, masterLayer);
    ras.push(hull.ra0);
    decs.push(hull.dec0);
  }

  // Pick the "middle" point if there are multiple master hulls (e.g. QA).
  // This fails if the hulls straddle ra=0/360.
  //
  let ra0 = 0;
  let dec0 = 0;
  if (ras.length > 1) {
    ra0 = 0.5 + (Math.min(...ras) + Math.max(...ras));
    dec0 = 0.5 + (Math.min(...decs) + Math.max(...decs));
  } else {
    ra0 = ras[0];
    dec0 = decs[0];
  }
  goToRaDec({ra: ra0, dec: dec0}, display);
}

// The following is OLD
//
// Returns information useful when the display is being created,
// but not when the regions are being re-created.
//
// TODO: rework this now moving knowledge into JS
//
function addRegionToJS9(img, ensemble, masterid, stack, winid,
                        regions) {

  /* can have multiple components for a stack */
  let linestyle;
  let nstackhulls = 0;

  const js9win = {display: img};

  const hulls = [];
  for (let i = 0; i < regions.stacks.length; i++) {
    let shull = regions.stacks[i];
    if (shull.stack === stack) {

      hulls.push(shull);

      /* experiment with making stack-level hulls "fixed", but
         I can imagine this could be annoying at times (useful
         at others) */
      linestyle = [1];  /* how to change this */

      add_hull_to_js9(shull,
                      {color: 'orange',
                       strokeDashArray: linestyle,
                       changeable: false,
                       tags: 'stack'},
                      js9win);

      nstackhulls += 1;
    }
  }

  // Ensure that the store is updated
  if (hulls.length > 0) {
    regionstore[stack] = hulls;
  } else {
    delete regionstore[stack];
  }

  /*
   * Add master hull last (so it appears on top).
   * Trying to add the right level of edit-ability to the polygon.
   *
   * Note that points can be moved and deleted even with these options
   * set.
   */

  if (typeof regions.master !== "undefined") {

    add_hull_to_js9(regions.master,
                    {movable: false,
                     rotatable: false,
                     resizable: false,  // do we want this true?
                     tags: 'master'},
                    js9win);

    regionstore.master = [regions.master];
  } else {
    delete regionstore.master;
  }

  const out = {nhulls: nstackhulls};
  if (typeof(regions.center) !== "undefined") {
    out.center = regions.center;
  }
  return out;
}

// Add the number of stacks to the JS9 window.
function updateJS9StackCounter(winid, n) {
  document.getElementById(winid + 'stackid')
    .innerHTML += "<span style='float: right;'>Stack hulls: " +
                  n.toString() + "</span>";
}

/*
 * HTML code for the JS9 display.
 *
 * Note that the event listeners are added in the
 * setupJS9 callback.
 */
function js9_display_html(stack, stacknum, id) {
  let html = "<div class='stackid' id='" + id + "stackid'>" +
    "<span style='float: left;'>Stack: " + stacknum + ": " +
    stack + "</span>" +
    "</div>" +
    "<div class='useropts'>";

  let name = id + "sigma";
  html += "<div class='blur'><span>Blur</span>";
  const sigmas = [0, 1, 2, 3, 4];
  for (let i = 0; i < sigmas.length; i++) {
    /* need to have label after the input for the CSS */
    const sigma = sigmas[i].toString();
    var l = id + 'sigma' + sigma;
    html += "<input class='sigma' name='" + name + "' type='radio' ";
    html += "id='" + l + "' value='" + sigma + "'"
    if (i === 0) { html += " checked"; }
    html += ">";
    html += "<label for='" + l + "'>" + sigma + "</label>";
  }
  html += "</div>";

  html += "<div class='reload'>";
  html += "<button id='" + id;
  html += "ReloadRegions'>Reload Regions</button>";
  html += "</div>";

  const all_bins = [1, 2, 4, 8, 16, 32, 64, 128];
  let def_binsize, def_start;
  if (stack.startsWith('hrc')) {
    def_binsize = 64;
    def_start = 1;
  } else {
    def_binsize = 8;
    def_start = 0;
  }
  const def_bins = all_bins.slice(def_start, def_start + 7);

  name = id + "binsize";
  html += "<div class='rebin'><span>Bin</span>";
  for (let i = 0; i < def_bins.length; i++) {
    const binsize = def_bins[i].toString();
    const l = id + 'binsize' + binsize;
    html += "<input class='binsize' name='" + name + "' type='radio' ";
    html += "id='" + l + "' value='" + binsize + "'";
    if (def_bins[i] === def_binsize) {
      html += " checked";
    }

    html += ">";
    html += "<label for='" + l + "'>" + binsize + "</label>";
  }
  html += "</div>";

  html += "<div class='zoom'>Zoom: ";
  html += "<button class='zoomin'>In</button>";
  html += "<button class='zoomout'>Out</button>";
  html += "</div>";

  html += "<div class='showable'>Toggle: ";
  html += "<button id='" + id + "ShowPanner'>Panner</button>";
  html += "</div>";

  html += "</div>";

  html += "<div class='JS9Menubar' id='" + id + "Menubar'></div>";
  html += "<div class='JS9' id='" + id + "'></div>'";
  html += "<div class='JS9Colorbar' id='" + id + "'></div>'";

  return html;
}


// "unique" id for new JS9 windows
//
//
var idctr = 1;
function js9_id(stack) {
  // return stack;

  var retval = idctr.toString();
  idctr += 1;
  return "js9win" + retval;
}

// Applies an energy filter for ACIS data to try and match the hull(s).
//
function showInJS9(stackval) {
  if (stackval.trim() === '') { return; }

  const toks = stackval.split(',');
  if (toks.length != 2) {
    alert("Internal error: stackval=[" + stackval + "]");
    return;
  }

  const stack = toks[0];
  const stacknum = toks[1];

  const winid = js9_id(stack);
  const opts = {onload: finalizeJS9Display('"""

    out += ensemble + "', " + str(masterid) + ", stack, winid)}"
    out += """
  if (stack.startsWith('acis')) {
    opts.bin = 8;
    opts.xdim = 8192;
    opts.ydim = 8192;
    // does this actually do anything?
    opts.filter = enbands[stack];
  } else {
    opts.bin = 64;
    opts.xdim = 16384;
    opts.ydim = 16384;
  }
  opts.id = winid;
  return JS9.LoadWindow('/evt3/' + stack, opts, 'light',
                        js9_display_html(stack, stacknum, opts.id));
}

/*** callback code for handling region changes ***/

/*
 * Let's see how well it works with the onchange handler.
 * If not, could we chain on the mouseup handler?
 */

// The code adds tags for master and stack which could be used
// to distinguish regions, although stack polygons should not
// be changeable.
//
// What about regions that the user has added directly?
//
function handleRegionChange(img, action) {
  console.log(">> mode=" + action.mode + " tags=" + action.tags);

  // For now, we only care about master tags;
  // this may need to be tweaked if we use the tag name as an id
  //
  if (action.tags[0] !== "master") { return; }

  if (action.mode === "update") {
    broadcastMasterUpdate(img, action);
  } else if (action.mode === "remove") {
    broadcastMasterDelete(img, action);
  }
}

// Tell the other windows about the new polygon.
//
function broadcastMasterUpdate(img, action) {

  // This is the window that the user has changed, so we don't want
  // to change this window.
  const js9win = img.display.id;
  const baseWin = {display: js9win};

  // Convert the polygon into a convex hull (if necessary).
  //
  let chull_sky = window.convexHull(action.pts);
  let chull_eqpos = [];
  for (let sky of chull_sky) {
    chull_eqpos.push(JS9.PixToWCS(sky.x, sky.y, baseWin));
  }

  // Now that the convex-hull version appears in a different layer,
  // the name can match the master name.
  //
  // const convexName = 'convex';
  const convexName = 'master';

  const convexOpts = {color: 'cyan',
                      strokeDashArray: [5, 3],
                      changeable: false,
                      evented: false,
                      tags: convexName};

  // Since using LightWindows, can look for div.dhtmlwindow
  // containers, and the knowledge that the id for this is
  // 'd' + id-of-js9-div
  //
  // Although this doesn't work as soon as the user creates a panner
  // or other related window. These extra divs appear to have
  // ids like js9win2_JS9Panner_lightDiv
  // as opposed to djs9win2, which we want. So for now just
  // reject anything with an underscore.
  //
  // The convex hull is
  //   a) drawn in all windows, including the one being edited
  //   b) drawn first, so that it hopefully appears below the
  //      master hull, so that the master hull can still be
  //      selected/edited.
  //

  const divs = document.getElementsByClassName('dhtmlwindow');
  for (let div of divs) {
    if (div.id.includes('_')) {
      continue;
    }
    const owin = div.id.substring(1);

    const imname = {display: owin};
    let hdl = JS9.GetImage(imname);

    // stop propogating the onchange signal
    hdl.params.xeqonchange = !hdl.params.xeqonchange;

    // Do we have a convex hull already?
    // let rs = JS9.GetRegions(convexName, imname);
    let rs = JS9.GetShapes(convexLayer, convexName, imname);
    if (rs.length === 0) {
      var ras = [];
      var decs = [];
      for (let wcs of chull_eqpos) {
        ras.push(wcs.ra);
        decs.push(wcs.dec);
      }
      add_hull_to_js9({ra: ras, dec: decs}, convexOpts,
                      imname, convexLayer);
    } else {
      const hullpts = [];
      for (let wcs of chull_eqpos) {
        hullpts.push(JS9.WCSToPix(wcs.ra, wcs.dec, imname));
      }

      // JS9.ChangeRegions(convexName, {pts: hullpts}, imname);
      JS9.ChangeShapes(convexLayer, convexName, {pts: hullpts},
                       imname);
    }

    // Only adjust the master polygon if this is not the
    // window the user is changing.
    if (owin !== js9win) {

      // Need to convert to image coordinates
      const pts = [];
      for (let wcs of action.wcspts) {
        pts.push(JS9.WCSToPix(wcs.ra, wcs.dec, imname));
      }

      // JS9.ChangeRegions('master', {pts: pts}, imname);
      JS9.ChangeShapes(masterLayer, 'master', {pts: pts}, imname);
    }

    hdl.params.xeqonchange = !hdl.params.xeqonchange;

    // Ensure the regions layer is where the user action happens
    JS9.ActiveShapeLayer(masterLayer, imname);
  }

}

// TODO: this copies from broadcastMasterUpdate, so need
//       to refactor
//
// Deletes the master shape from both the region and
// convex layers. It does not delete the original version.
//
function broadcastMasterDelete(img, action) {

  const js9win = img.display.id;

  // Since using LightWindows, can look for div.dhtmlwindow
  // containers, and the knowledge that the id for this is
  // 'd' + id-of-js9-div
  //
  // Although this doesn't work as soon as the user creates a panner
  // or other related window. These extra divs appear to have
  // ids like js9win2_JS9Panner_lightDiv
  // as ooposed to djs9win2, which we want. So for now just
  // reject anything with an underscore.
  //
  const divs = document.getElementsByClassName('dhtmlwindow');
  for (let div of divs) {
    if (div.id.includes('_')) {
      continue;
    }
    const owin = div.id.substring(1);
    if (owin === js9win) {
      continue;
    }

    const imname = {display: owin};
    let hdl = JS9.GetImage(imname);

    // stop propogating the onchange signal
    hdl.params.xeqonchange = !hdl.params.xeqonchange;

    // JS9.RemoveRegions('master', imname);
    for (let layer of [convexLayer, masterLayer]) {
      JS9.RemoveShapes(layer, 'master', imname);
    }

    hdl.params.xeqonchange = !hdl.params.xeqonchange;
  }

  // Delete the convex-hull version in the original window.
  JS9.RemoveShapes(convexLayer, 'master', {display: js9win});

}

if (JS9.Regions.opts.onchange !== null) {
  alert("Overwriting onchange handler!");
}

JS9.Regions.opts.onchange = handleRegionChange;

/* play with websamp */

var sampConnection;

// UGLY
var sampBaseUrl = window.location.href.toString()
                     .replace(new RegExp("/[^/]*/[^/]*/[^/]*$"), "");

function showInSAMP(stack) {
  if (stack.trim() === '') { return; }
  const url = sampBaseUrl + "/evt3/" + stack;
  console.log("*** SAMP : " + url);

  sampConnection.runWithConnection((connection) => {
    var msg = new samp.Message("image.load.fits",
                               {"url": url, "name": stack});
    connection.notifyAll([msg]);
  });
}

function sampIsAvailable(flag) {
  document.getElementById("sampcontrols").hidden = !flag;
}

function finalize() {
  if (typeof sampConnection !== "undefined") {
    sampConnection.unregister();
    sampConnection = undefined;
  }
}

function initialize() {

  /*** For the moment comment out the samp code, as the
       continual checking is making some unrelated
       tasks hard to debug/work on.
  if (typeof sampConnection === "undefined") {
    sampConnection = new samp.Connector("CHS Sender");
    sampConnection.onHubAvailability(sampIsAvailable, 2000);
  }
  ***/
  sampIsAvailable(false);


  document.getElementById('imgscale').
    addEventListener("change", (e) => { setScaling(e.target.value); });

  document.getElementById('showinsamp').
    addEventListener("change", (e) => { showInSAMP(e.target.value); });

  document.getElementById('showinjs9').
    addEventListener("change", (e) => { showInJS9(e.target.value); });

  // Set up the button handlers, if we have any
  //
  if (npages > 1) {
    for (let i = 1; i <= npages; i++) {
      document.getElementById("page" + i.toString())
        .addEventListener("change", (e) => { setPage(e.target.value); });
    }
  }

  setPage(1);
}
</script></head>"""

    # out += "<body onload='initialize();' onunload='finalize();'>"
    out += "<body onload='initialize();'>"

    out += "<div id='infobar'>"

    out += "<div id='home'>"
    out += "<a href='/'>Ensemble list</a>"
    out += "<br>"
    out += "<a href='/{0}'>{0}</a>".format(ensemble)
    out += "<br>"

    # button list to get to the other master hulls
    for h in info['masters']:
        mval = int(h['masterid'])
        if h['masterid'] == mid:
            out += "<span class='hulllink'>{}</span>".format(mval)
            continue

        out += "<a class='hulllink' href='/{}/{}/{}'>{}</a>".format(ensemble,
                                                                    revstr,
                                                                    h['masterid'],
                                                                    mval)

    out += "</div>"

    out += "<div class='info'>"
    out += "<span class='label'>{}</span>".format(ensemble)

    out += "<span>Hull:</span>"
    out += "<span id='hullid'>{}</span>".format(mid)

    out += "<span>Number components:</span>"
    out += "<span id='ncomponents'>{}</span>".format(hull['ncpts'])

    out += "<span>Version:</span>"
    out += "<span>{}</span>".format(revstr)
    out += "</div>"

    # TODO: is select the best UI option for the 'select a stack to
    #       display' option?

    # Try WebSamp
    #
    out += "<div id='sampcontrols'>SAMP: "
    out += "<select id='showinsamp'>"
    out += "<option value='' selected></option>"
    for stk in ordered_stks:
        out += "<option value='{0}'>{1:03d}: {0}</option>".format(stk,
                                                                  info['stackmap'][stk])
    out += "</select>"
    out += "</div>"

    # Add in a JS9 link to select the given stack.
    #
    out += "<div id='js9controls'>"
    out += "JS9: <select id='showinjs9'>"
    out += "<option value='' selected></option>"
    for stk in ordered_stks:
        out += "<option "
        out += "value='{0},{1:03d}'>{1:03d}: {0}".format(stk,
                                                         info['stackmap'][stk])
        out += "</option>"

    out += "</select>"
    out += "</div>"

    out += "<div id='summary'><select id='imgscale' class='button'>"
    out += "<option value='log10' selected>log 10</option>"
    out += "<option value='sqrt'>square root</option>"
    out += "<option value='none'>linear</option>"
    out += "</select>"

    # list the pages, if there's more than one.
    #
    # the pages could just be a list of buttons, in a grid, which
    # change to indicate the selected one, but that's more work.
    #
    if hull['npages'] > 1:
        out += "<table><thead><tr><th>Page</th><th>Select</th>"
        out += "</tr></thead><tbody>"

        for p in range(1, hull['npages'] + 1):
            pageid = 'page{}'.format(p)
            out += "<tr><td><label for='{}'>{}<label></td>".format(pageid, p)
            out += "<td><input name='pagechoice' type='radio' id='{}' ".format(pageid)
            out += "value='{}'".format(p)
            if p == 1:
                out += "checked"
            out += "></td></tr>"

        out += "</tbody></table>"

    out += "</div>"

    # notes
    out += "<div id='notes'><p class='label'>Notes</p>"
    out += "<textarea cols=32 rows=10 id='usercontent'></textarea>"
    out += "<br><button id='saveusercontent' class='button'>"
    out += "Save notes</button></div>"

    # user action
    #
    # assumes the useraction field matches these options.
    #
    out += "<div id='action'><p>User action:</p>"
    if is_latest:
        out += "<select required='true' id='choice' class='button'>"
        for n, v in [("", ""), ("Accept", "accept"),
                     ("Delete Master", "delete"),
                     ("Manual", "manual")]:
            out += "<option value='{}'"
            if v == hull['useraction']:
                out += " selected"

            out += ">{}</option>".format(v, n)

        out += "</select>"
    else:
        try:
            action = {'': 'ERROR: no action recorded',
                      'accept': 'Accept',
                      'delete': 'Delete Master',
                      'manual': 'Manual'}[hull['useraction']]
        except KeyError:
            print("LOG: useraction unknown: {}".format(hull['useraction']))
            action = 'ERROR: unrecognized action ' + \
                '<{}>'.format(hull['useraction'])

        out += "<p>{}</p>".format(action)
    out += "</div>"

    # overview image
    out += "</div><div id='overviewbar'>"
    out += "<div id='username'>{}</div>".format(os.getlogin())

    out += "<img id='pageview'>"
    out += "</div>"

    out += "</body></html>"
    return 200, out


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
        print("LOG: error reading [{}]: {}".format(infile, exc))
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
        upath = urlparse.urlparse(self.path)
        path = upath.path

        # Strip the leading /
        if path != "":
            path = path[1:]

        if path in ["", "index.html"]:
            status, cts = create_index_page(datadir)
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
                status, cts = create_ensemble_index_page(datadir,
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
                status, cts = create_master_hull_page(datadir,
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
            print("LOG: invalid master id: [{}]".format(toks[2]))
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
            print("LOG: unknown JS9 suffix [{}]".format(infile))
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
            print("LOG: failed glob={}".format(pat))
            self.send_error(404)
            return

        # Should pick the highest version; for now pick the first
        #
        infile = matches[0]

        if use_gzip:
            try:
                cts = open(infile, 'r').read()
            except IOError as exc:
                print("LOG: error reading [{}]: ".format(infile) +
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


def serve(userdir, webdir, datadir, port=8070,
          evt3dir='/data/L3/chs_master_match/input/stkevt3'):

    for dirname in [userdir, webdir, datadir, evt3dir]:
        if not os.path.isdir(dirname):
            raise IOError("Not a directory: {}".format(dirname))

    store = {'userdir': userdir,
             'webdir': webdir,
             'datadir': datadir,
             'evt3dir': evt3dir}

    server_address = ('', port)
    httpd = CHSServer(store, server_address, CHSHandler)
    print("Starting server on http://localhost:{}/".format(port))
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

    # Store the user's files in the current working directory
    #
    userdir = os.getcwd()

    serve(userdir=userdir, webdir=webassets,
          datadir=datadir, port=port)
