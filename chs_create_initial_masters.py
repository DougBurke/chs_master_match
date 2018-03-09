#!/usr/bin/env python

"""
Usage:

  ./chs_create_initial_masters.py ensemblefile ensemble outdir
      --mrgsrc3dir <dirname>
      --compzero n

Aim:

Create the initial merged-hull table for the Convex-Hull Master match
of an ensemble. The output directory (outdir) must not exist when the
script is run.

The ensemblefilefile argument is a file containing the mapping from
ensemble to stack. It should be readable by crates and contain the
columns ensemble and stack. It is expected to be the file created by
get_ensemble_to_stack_mapping.py. An alternative approach would be to
read in the ensemble mst3 file, which would give the ensemble name and
the stack ids.

If there are no hulls in the ensemble then the file outdir/NOHULLS
is created (it is an error if this file already exists) and the script
exits with a status of 0.

The compzero value is added to the component numbers of the stack-level
hulls to create the COMPONENT column of the HULLMATCH block. It is
added as the COMPZERO keyword to the header, and defaults to 9000 when
not given. It refers to the COMPONENT value of a stack-level hull
with a zero component value (they are 1-based).

If there are any hulls then the following file is created in the
directory outdir:

  master_hulls.<ensemble>.v001.fits

This file is created even if all the master hulls are in a QA state.

If there are any qa cases then they are written to

  qa.<ctr>.v001.fits

where <cpt> is the Master_Id value (written out as 3-digit 0-padded
integer) of the master source: there will be multiple polygons in the
file (as this is the only current send-to-QA case).

There are also DS9 ascii region files (in celestial coordinates)
created for each hull: stack, master, and qa.

  master.<cpt>.v001.reg   - unless there is no master hull
  qa.<cpt>.v001.reg
  stack.<stack>.<stackcpt>.v001.reg

TODO: there should be some visual distinction between these regions
      when displayed in ds9.

Notes
-----

The Master_Id values are *NOT* guaranteed to be consistent from run to
run, even if the data has not changed.

It does not deal nicely with merged hulls that are not 'okay' or 'qa';
e.g. if there is an error creating the hull.

"""

import itertools
import os
import time

import numpy as np

import pycrates
import region

import chs_utils as utils
import chs_identify_master_hulls as identify
import chs_merge_hulls as merge

NODATAFILE = "NOHULLS"

help_str = "Create the review products for CHS in this ensemble."


def find_ensemble_stacks(ensemblefile, ensemble):
    """Return the stacks associated with the ensmeble.

    Parameters
    ----------
    emsemble : str
        The ensemble name.
    ensemblefile : str
        The file should contain columns ensemble and stack, and is
        used to find what stacks to look for.

    Returns
    -------
    stacks : list of str
        The stack names.

    Raises
    ------
    IOError
        If the ensemble is not known about
    """

    # I could use an 'ensemble=<...>' filter, but there are known
    # DM string-filtering bugs in CIAO 4.9 which mean that we
    # have to do the filtering manually.
    #
    infile = "{}[cols ensemble, stack]".format(ensemblefile)
    cr = pycrates.read_file(infile)

    # do not assume the input file is sorted by ensemble
    out = []
    for ens, stack in zip(cr.ensemble.values,
                          cr.stack.values.copy()):
        if ens != ensemble:
            continue

        out.append(stack)

    if out == []:
        raise IOError("No stacks found for " +
                      "ensemble {} in {}".format(ensemble,
                                                 ensemblefile))

    return out


def read_hulls(stack, mrgsrc3dir):
    """Read in the hulls for the stack.

    This only processes hulls with a STATUS of 0.

    Parameters
    ----------
    stack : str
        The stack name.
    mrgsrc3dir : str
        The location of the mrgsrc3 files.

    Returns
    -------
    hulls : list of dict
        The hull information. Each element is a hull, containing
        keys: stack, component, transform, pos, eqpos, area, infile.
        The coordinates are 2D (2, npts) and have been filtered
        to ensure that there are no non-finite numbers, and that
        the region is closed. The area is the area of the polygon
        in square arcsec.

    Raises
    ------
    IOError
        If the STACK_ID keyword in the MEXTSRC block of the mrgsrc3
        files does not match the stack argument.
    """

    mrgsrc3file = utils.find_mrgsrc3(stack, mrgsrc3dir)
    infile = "{}[MEXTSRC][status=0]".format(mrgsrc3file)

    cr = pycrates.read_file(infile)
    stack_key = cr.get_key_value('STACK_ID')
    if stack != stack_key:
        raise IOError("STACK_ID mismatch: " +
                      "{} vs {} in\n{}".format(stack, stack_key,
                                               infile))

    # doesn't really save much work
    if cr.get_nrows() == 0:
        return []

    transform = cr.get_transform('EQSRC').copy()

    # One pixel in arcsec
    pixsize = 3600.0 * transform.get_parameter_value('CDELT')[1]

    # For now it is easier to deal with the MAN_CODE as a byte
    # rather than a bit array, so switch.
    #
    cr.MAN_CODE.convert_bits_to_bytes()

    out = []
    for cpt, eband, mancode, pos, eqsrc in zip(cr.COMPONENT.values.copy(),
                                               cr.EBAND.values.copy(),
                                               cr.MAN_CODE.values.copy(),
                                               cr.POS.values.copy(),
                                               cr.get_column('EQSRC').values.copy()):

        pos = utils.validate_polygon(pos, report=True,
                                     label="{} {}".format(stack,
                                                          cpt))
        eqsrc = utils.validate_polygon(eqsrc, report=False)
        assert pos.shape == eqsrc.shape

        cstr = ",".join(["{}d".format(p)
                         for p in pos.T.flatten()])
        regstr = "polygon({})".format(cstr)
        reg = region.regParse(regstr)
        area_pixels = region.regArea(reg)

        out.append({'stack': stack,
                    'component': cpt,
                    'infile': infile,
                    'transform': transform,
                    'area': area_pixels * pixsize * pixsize,
                    'eband': eband,
                    # although converted to a byte, it is still a vector
                    # column, although only of length 1
                    'man_code': mancode[0],
                    'pos': pos,
                    'eqpos': eqsrc})

    return out


def no_hulls(outdir):
    """This ensemble has no hulls.

    Create the output directory and add in the marker file.

    Parameters
    ----------
    outdir : str
        The output directory, which will be created by the routine.
    """

    # See
    # https://stackoverflow.com/questions/1158076/implement-touch-using-python
    # which is a bit OTT for our needs.
    #
    os.mkdir(outdir)
    outfile = os.path.join(outdir, NODATAFILE)
    open(outfile, 'a').close()


# Should send in the time value
_header = {}


def add_standard_header(cr, creator=None, revision=1):

    if 'timestr' not in _header:
        _header['timestr'] = time.strftime("%Y-%m-%dT%H:%M:%S")

    timestr = _header['timestr']

    hdrvals = []
    if creator is not None:
        hdrvals.append(('CREATOR', creator,
                        'tool that created this output'))

    hdrvals.extend([('DATE', timestr,
                     'Date and time of file creation'),
                    ('CHSVER', revision,
                     'The revision number for the file')
                    ])

    add_header(cr, hdrvals)


def add_header(cr, hdrvals):
    """Add the header keywords to the crate.

    Parameters
    ----------
    cr : pycrates.Crate
    hdrvals : sequence of (name, value, desc) triples
    """

    for name, value, desc in hdrvals:
        key = pycrates.CrateKey()
        key.name = name
        key.value = value
        key.desc = desc
        cr.add_key(key)


def add_col(cr, name, values,
            desc=None,
            unit=None):
    """Add a column to the crate.

    Parameters
    ----------
    cr : pycrates.TABLECrate
    name : str
    values : array_like
    desc : str or None, optional
        The description.
    unit : str or None, optional
        The units for the column.

    """

    col = pycrates.CrateData()
    col.name = name
    col.values = values
    if desc is not None:
        col.desc = desc
    if unit is not None:
        col.unit = unit

    cr.add_column(col)


def write_hulls(ensemble, outfile, hullcpts, hullareas, outlines,
                hull_store,
                stacks=None,
                compzero=9000,
                revision=1,
                creator=None):
    """Create the "CHS mst3" file.

    Parameters
    ----------
    ensemble : string
        The ensemble value, written to the header as the ENSEMBLE
        keyword.
    outfile : string
        This file is overwritten if it exists.
    hullcpts : sequence of (stack, cpt) sequences
        The stack-level hulls that form each master hull.
    hullareas : dict
        The stack hull areas, in arcsec^2, where the key is the
        pair (stack, component) and the value is the area.
    outlines : list of dict
        The corresponding outline (can be an error state).
    hull_store : dict
        The keys are (stack, component) and the values are the
        stack-level data read from the mrgsrc3 file.
    stacks : list of str or None, optional
        The stacks that form this ensemble (whether or not they
        have a stack-level convex hull). If given then the
        output file will contain STKIDNUM and STKIDxxx values.
    compzero : int, optional
        The value of the COMPONENT column in the HULLMATCH block
        for a stack-level hull which has a component value of 0
        (from the mrgsrc3 MEXTSRC block); note that there are no
        such component values since they are 1 based. Must be >= 0.
    revision : int
        The value to write out to the header as the CHSVER
        keyword.
    creator : None or str, optional
        The name to use for the CREATOR field in the header.
    """

    assert len(hullcpts) == len(outlines), \
        'len = {} vs {}'.format(len(hullcpts), len(outlines))

    extra_hdr = [('ENSEMBLE', ensemble, 'The ensemble'),
                 ('COMPZERO', compzero, 'The COMPONENT value for cpt=0')]
    if stacks is not None:
        stacks = sorted(list(stacks))
        for i, stack in enumerate(stacks):
            extra_hdr.append(('STKID{:03d}'.format(i),
                              stack,
                              'Member of the ensemble'))

        extra_hdr.append(('STKIDNUM', len(stacks),
                          'Number of stacks in ensemble'))

    ds = pycrates.CrateDataset()

    cr = pycrates.TABLECrate()
    cr.name = 'HULLMATCH'

    add_standard_header(cr, creator=creator, revision=revision)
    add_header(cr, extra_hdr)

    # What data do we want to store?
    mid = []
    nhulls = []
    stks = []
    cpts = []
    mtypes = []
    areas = []
    ebands = []
    mancodes = []
    for i, hcpts in enumerate(hullcpts):

        m = i + 1
        nh = len(hcpts)

        # sort to try and make the file easier to scan by the user
        for key in sorted(hcpts):
            assert key in hull_store

            stack, cpt = key

            mid.append(m)
            nhulls.append(nh)
            stks.append(stack)
            cpts.append(cpt + compzero)
            mtypes.append('Unambiguous')
            areas.append(hullareas[key])
            ebands.append(hull_store[key]['eband'])

            # would like to keep as a bit array
            mancodes.append(hull_store[key]['man_code'])

    # NOTE: the NHULLS column might be better in the next block,
    #       but it can be useful to know how many stack-level
    #       hulls there are in a master when looking at this data.
    #
    add_col(cr, 'Master_Id', mid,
            desc='This is an internal number, do not expose')
    add_col(cr, 'NHULLS', nhulls,
            desc='The number of stack-level hulls in the master')
    add_col(cr, 'STACKID', stks)
    add_col(cr, 'COMPONENT', cpts,
            desc='Offset by COMPZERO from MEXTSRC component value')
    add_col(cr, 'Match_Type', mtypes)
    add_col(cr, 'AREA', areas,
            unit='arcsec**2',
            desc='Area of hull excluding pixel-mask filter')
    add_col(cr, 'EBAND', ebands,
            desc='Energy band of hull')
    add_col(cr, 'MAN_CODE', mancodes,
            desc='Copied from MEXTSRC block (converted to int)')

    ds.add_crate(cr)

    cr = pycrates.TABLECrate()
    cr.name = 'HULLLIST'

    add_standard_header(cr, creator=creator, revision=revision)
    add_header(cr, extra_hdr)

    # What is the maximum number of points in a hull?
    #
    nmax = 0
    for outline in outlines:
        # should we skip the QA case?
        if outline['status'] != 'okay':
            continue

        eqpos = outline['eqpos']
        nmax = max(nmax, eqpos.shape[1])

    # It is okay for nmax=0 - if all master hulls are QA cases - so
    # only want to warn if we have an unusually-small polygon.
    #
    if nmax < 3:
        if nmax > 0:
            print("WARNING: max number of vertices in a " +
                  "hull={}".format(nmax))
        nmax = 3

    # Some of these columns are a single value, but put in here
    # so we can see the proposed file structure, and to make
    # downstream processing a bit easier (i.e. not having to worry
    # about whether a file exists).
    #
    mid = []
    status = []
    base_stack = []
    man_reg = []
    man_match = []

    nhulls = len(outlines)
    eqpos = np.full((nhulls, 2, nmax), np.nan,
                    dtype=np.float64)

    # Ugh: loop getting messy
    nstkhull = [len(c) for c in hullcpts]
    nvertex = []

    for i, outline in enumerate(outlines):

        # TODO: could calculate the number of stacks that contribute/
        #       number with no data

        mid.append(i + 1)
        status.append(outline['status'])

        man_reg.append(False)
        man_match.append(False)

        bs = outline['base_stack']
        if bs is None:
            bs = "NONE"

        base_stack.append(bs)

        if outline['status'] == 'okay':
            vs = outline['eqpos']
            assert vs.ndim == 2
            assert vs.shape[0] == 2
            npts = vs.shape[1]

            nvertex.append(npts)
            eqpos[i, :, :npts] = vs

        else:
            nvertex.append(0)

            print("Skipping data for " +
                  "Match_Id={} status={}".format(i + 1,
                                                 outline['status']))

    add_col(cr, 'Master_Id', mid,
            desc='This is an internal number, do not expose')
    add_col(cr, 'STATUS', status,
            desc='Did the master-match work?')
    add_col(cr, 'BASE_STK', base_stack,
            desc='The stack used for SKY coord system, or NONE')

    add_col(cr, 'MANMATCH', man_match,
            desc='Has the selection of stack-level hulls been changed')
    add_col(cr, 'MANREG', man_reg,
            desc='Has the region been changed manually')

    add_col(cr, 'NVERTEX', nvertex,
            desc='The number of vertexes in the closed hull')

    add_col(cr, 'NSTKHULL', nstkhull,
            desc='The number of stack hulls that were combined')

    # NOTE: there is no POS coordinate column in this block
    col = pycrates.create_vector_column('EQPOS', ['RA', 'DEC'])
    col.desc = 'The master hull vertices'
    col.unit = 'degree'
    col.values = eqpos
    cr.add_column(col)

    ds.add_crate(cr)

    ds.write(outfile, clobber=True)
    print("Created: {}".format(outfile))


def dump_qa(ensemble, outdir, ctr, outline,
            creator=None, revision=1, color='green'):
    """Dump the QA data to a region-like FITS file

    The output files are qa.<ctr>.v<version>.[fits|reg], where the
    integers are written as 3-character, zero-padded values.

    Parameters
    ----------
    ensemble : string
        The ensemble value, written to the header as the ENSEMBLE
        keyword.
    outdir : string
        This output directory
    ctr : int
        Used for the file name and added to the header as the
        HULLCPT keyword.
    outline : dict
        The hull data to write out.
    creator : None or str, optional
        The name to use for the CREATOR field in the header.
    revision : int
        The value to write out to the header as the CHSVER
        keyword.
    color : str
        The color for the DS9 region files.

    Notes
    ------
    I'd like to add a transform to POS but it is a little-bit involved
    and I don't have time to dig up the code, so just add an explicit
    column instead of a virtual one. It should be the case that the
    base stack is the same for each row, but leave as a column rather
    than move to a header keyword for now.

    """

    assert outline['status'] == 'qa'
    assert outline['eqpos'] is not None
    assert outline['pos'] is not None
    assert outline['base_stack'] is not None

    outfile = os.path.join(outdir,
                           'qa.{:03d}.v{:03d}.fits'.format(ctr,
                                                           revision))

    cr = pycrates.TABLECrate()
    cr.name = 'QACASE'

    add_standard_header(cr, creator=creator, revision=revision)
    add_header(cr, [('ENSEMBLE', ensemble,
                     'The ensemble'),
                    ('HULLCPT', ctr,
                     'The Master_Id of the hull')])

    eqpos = outline['eqpos']
    pos = outline['pos']
    assert len(pos) == len(eqpos)

    ncpts = len(pos)

    add_col(cr, 'COMPONENT', np.arange(1, ncpts + 1))
    add_col(cr, 'SHAPE', ['Polygon'] * ncpts)
    add_col(cr, 'BASE_STK', [outline['base_stack']] * ncpts,
            desc='The stack used for SKY coord system, or NONE')

    nvertex = []
    for poly in pos:
        x = poly[0]
        y = poly[1]
        xidx = np.isfinite(x)
        yidx = np.isfinite(y)
        assert (xidx == yidx).all()
        nvertex.append(xidx.sum())

    add_col(cr, 'NVERTEX', nvertex,
            desc='The number of vertexes in the closed hull')

    nmax = max(nvertex)
    pos_out = np.full((ncpts, 2, nmax), np.nan, dtype=np.float64)
    eqpos_out = np.full((ncpts, 2, nmax), np.nan, dtype=np.float64)

    for i, dvals in enumerate(zip(pos, eqpos)):

        npts = nvertex[i]
        phys, cel = dvals
        pos_out[i, :, :npts] = phys
        eqpos_out[i, :, :npts] = cel

    col = pycrates.create_vector_column('POS', ['X', 'Y'])
    col.desc = 'The master hull vertices'
    col.unit = 'pixel'
    col.values = pos_out
    cr.add_column(col)

    col = pycrates.create_vector_column('EQPOS', ['RA', 'DEC'])
    col.desc = 'The master hull vertices'
    col.unit = 'degree'
    col.values = eqpos_out
    cr.add_column(col)

    cr.write(outfile)
    print("Created: {}".format(outfile))

    outfile = os.path.join(outdir,
                           'qa.{:03d}.v{:03d}.reg'.format(ctr,
                                                          revision))
    with open(outfile, 'w') as ofh:
        ds9_header(ofh, color=color)

        for i, cel in enumerate(eqpos):

            npts = nvertex[i]
            ostr = ds9_shape(cel[:, :npts])
            ostr += ' # text={{Id={} {}}}\n'.format(ctr,
                                                    i + 1)
            ofh.write(ostr)

    print("Created: {}".format(outfile))


def ds9_shape(cel):
    """Given celestial coordinates, return the polygon string.

    Parameters
    ----------
    cel : numpy array
        The shape must be (2, nvertex) and all vertices are output
        (i.e. it must already have been filtered for NaN).

    Returns
    -------
    poly : str
        The DS9 polygon representation.

    """

    if not np.all(np.isfinite(cel)):
        raise ValueError("Found non-finite value")

    celstr = ",".join(['{}d'.format(p) for p in cel.T.flatten()])
    return 'polygon({})'.format(celstr)


def ds9_header(ofh, color='green'):
    """Write out the starter for a DS9 region file."""

    ofh.write('# Region file format: DS9 version 4.1\n')
    ofh.write('global color={} '.format(color) +
              'dashlist=8 3 width=1 ' +
              'font="helvetica 10 normal roman" select=1 ' +
              'highlite=1 dash=0 fixed=0 edit=1 move=1 ' +
              'delete=1 include=1 source=1\n')
    ofh.write('fk5\n')


def write_stack_hull_as_ds9(hull, outdir, revision, color='green'):
    """Write the stack hull as a ds9 region file."""

    outfile = os.path.join(outdir,
                           'stack.{}.{}.v{:03d}.reg'.format(hull['stack'],
                                                            hull['component'],
                                                            revision))
    with open(outfile, 'w') as ofh:
        ds9_header(ofh, color=color)
        ostr = ds9_shape(hull['eqpos'])
        ostr += ' # text={{stack={} {}}}\n'.format(hull['stack'],
                                                   hull['component'])
        ofh.write(ostr)


# TODO: move to chs_merge_hulls?
#
def do_masters_overlap(m1, m2, transforms):
    """Do the two master hulls overlap?

    Parameters
    ----------
    m1, m2 : dict
        The two master hulls to check. The dicts are the return value
        from chs_merge_hulls.make_merged_hulls. The contents of these
        dicts may be changed by this routine.
    transforms : dict
        The keys are the stack and the values are the SKY -> EQ
        transform for that stack.

    Return
    ------
    overlap : bool
        Do the hulls overlap?

    Notes
    -----
    If a hull is in qa state then the check is made for each polygon
    (it would probably make sense to calculate the convex hull of the
    components and use that but that is  potential upgrade).

    Each outline looks like:

    {'status': 'okay', 'base_stack': 'acisfJ1711499m393614_001', 'eqpos': array([[ 258.35268129,  258.34706199,  258.31750321,  258.30620794,
         258.28499025,  258.18857561,  258.15412969,  258.33819263,
         258.35268129],
       [ -39.51102123,  -39.51978671,  -39.55924129,  -39.5691155 ,
         -39.58120352,  -39.5956472 ,  -39.45685734,  -39.45312316,
         -39.51102123]]), 'pos': array([[ 1868.5,  1900.5,  2068.5,  2132.5,  2252.5,  2796.5,  2988.5,
         1948.5,  1868.5],
       [ 4772.5,  4708.5,  4420.5,  4348.5,  4260.5,  4156.5,  5172.5,
         5196.5,  4772.5]])}

    """

    if m1['status'] == 'error' or m2['status'] == 'error':
        return False

    assert m1['pos'] is not None, 'master 1 POS is None'
    assert m2['pos'] is not None, 'master 2 POS is None'

    # We only care about stack for m1 and only in certain cases, but
    # check anyway (to validate assumptions).
    #
    stack1 = m1['base_stack']
    stack2 = m2['base_stack']
    assert stack1 in transforms, 'master 1 unknown transform'
    assert stack2 in transforms, 'master 2 unknown transform'

    # Need to convert to the same SKY coordinate system; it is not
    # clear to me whether hulls with the same base stack can overlap
    # so best check everything.
    #
    # For now in the situation where pos could be a 2D NumPy array
    # or a list of NumPy 2D arrays. It is not unreasonable to have
    # everything be a list here.
    #
    def listify(x):
        if type(x) == list:
            return x
        assert x.ndim == 2
        return [x]

    p1 = listify(m1['pos'])

    if stack1 == stack2:
        p2 = listify(m2['pos'])
    else:
        s2 = listify(m2['eqpos'])
        tr = transforms[stack1]
        p2 = []
        for s in s2:
            # Add in a dummy axis to simplify conversion, due to the
            # way the apply and invert methods work.
            #
            p = tr.invert(s[np.newaxis])
            p2.append(p[0])

    # Since we occasionally have to deal with multiple polygons per
    # master, ensure we do it for all cases.
    #
    def toqa(m):
        "Change status and pos/eqpos arrays"
        m['status'] = 'qa'
        m['pos'] = listify(m['pos'])
        m['eqpos'] = listify(m['eqpos'])

    # If any overlap we return True
    #
    h1s = [utils.make_region_string(h1) for h1 in p1]
    h2s = [utils.make_region_string(h2) for h2 in p2]
    for h1 in h1s:
        for h2 in h2s:
            if utils.check_overlap(h1, h2):
                toqa(m1)
                toqa(m2)
                return True

    return False


def process_ensemble(ensemblefile, ensemble, outdir,
                     mrgsrc3dir,
                     compzero=9000,
                     revision=1,
                     creator=None,
                     master_color='green',
                     stack_color='green',
                     qa_color='green'):
    """Process the ensemble to create the CHS review products.

    Parameters
    ----------
    ensemblefile : str
        The file should contain columns ensemble and stack, and is
        used to find what stacks to look for.
    ensemble : str
        The ensemble name.
    outdir : str
        The output directory, which will be created by the routine.
    mrgsrc3dir : str
        The directory name containing the mrgsrc3 files for the
        stacks. The names must match
        <stack>*mrgsrc3.fits* and there can only be one file per
        stack.
    compzero : int, optional
        The value of the COMPONENT column in the HULLMATCH block
        for a stack-level hull which has a component value of 0
        (from the mrgsrc3 MEXTSRC block); note that there are no
        such component values since they are 1 based. Must be >= 0.
    revision : int, optional
        The value to write out to the header as the CHSVER
        keyword.
    creator : str or None, optional
        If set it is used as the CREATOR keyword value in the output
        file.
    """

    if compzero < 0:
        raise ValueError("compzero must be >= 0, sent {}".format(compzero))

    if os.path.exists(outdir):
        raise IOError("The output directory already exists: {}".format(outdir))

    # validate input before creating the output
    #
    stacks = find_ensemble_stacks(ensemblefile, ensemble)
    nstacks = len(stacks)

    # I need to pass around eband,mancode info to write_hulls and
    # I am too lazy to re-architect the code, so I am adding a "global"
    # dict which contains the stack-level information indexed by
    # (STACKID, COMPONENT), since that is used in write_hulls.
    #
    hull_store = {}
    transform_store = {}

    hulls = []
    stacks_with_hulls = 0
    for stack in stacks:
        shulls = read_hulls(stack, mrgsrc3dir)
        if shulls == []:
            continue

        hulls.extend(shulls)
        stacks_with_hulls += 1

        for shull in shulls:
            key = (shull['stack'], shull['component'])
            assert key not in hull_store
            hull_store[key] = shull

        # Add in the stack transform
        transform_store[stack] = shulls[0]['transform']

    if hulls == []:
        print("No hulls were found in ensemble {}".format(ensemble))
        no_hulls(outdir)
        return

    os.mkdir(outdir)

    nhulls_stack = len(hulls)
    print("There are {} stack-level hulls ".format(nhulls_stack) +
          "in {} stacks".format(stacks_with_hulls))
    if stacks_with_hulls != nstacks:
        print("Total number of stacks: {}".format(nstacks))

    # Dump the stack-level hulls as a DS9 region file and extract
    # the hull areas.
    hull_areas = {}
    for hull in hulls:
        write_stack_hull_as_ds9(hull, outdir, revision=revision,
                                color=stack_color)

        hull_areas[hull['stack'], hull['component']] = hull['area']

    print("Created stack-level region files in {}".format(outdir))

    # What hulls overlap and what don't?
    #
    overlap_gr, singles = identify.find_overlap_graph(hulls)
    overlaps = identify.get_nodes(overlap_gr)

    master_hulls = []

    # Overlap hulls
    #
    seen = set([])
    for overlap in overlaps:
        master_hulls.append([(stack, cpt) for stack, cpt in overlap])
        for k in overlap:
            assert k not in seen, \
                'repeat occurence of {}'.format(k)
            seen.add(k)

    # Single hulls
    #
    for key in singles:
        master_hulls.append([key])

        assert key not in seen, \
            'repeat occurence of {}'.format(key)
        seen.add(key)

    nfound = len(seen)
    assert nhulls_stack == nfound, 'Lost (or gained) a hull'

    noverlap = len(overlaps)
    nsingle = len(singles)
    print("Found {} overlap and {} single hulls".format(noverlap,
                                                        nsingle))

    # Create the initial version of the master hulls and dump out
    # any QA cases here.
    #
    # The master hulls are created, then a check is made to
    # ensure they do not overlap (marking them as qa cases if
    # they do), then they are written out.
    #
    # THE OVERLAP CHECK IS NOT WRITTEN YET!
    #
    outlines = []
    for hullcpts in master_hulls:
        outline = merge.make_merged_hull(hullcpts, mrgsrc3dir)
        outlines.append(outline)

    noverlap = 0
    for m1, m2 in itertools.combinations(outlines, 2):
        if do_masters_overlap(m1, m2, transform_store):
            noverlap += 1

    if noverlap > 0:
        print("NOTE: found {} master overlaps -> QA".format(noverlap))

    for i, outline in enumerate(outlines):
        if outline['status'] != 'okay':
            print("status = {}".format(outline['status']))
        else:
            outfile = os.path.join(outdir,
                                   'master.{}.v{:03d}.reg'.format(i + 1,
                                                                  revision))
            with open(outfile, 'w') as ofh:
                ds9_header(ofh, color=master_color)
                ostr = ds9_shape(outline['eqpos'])
                ostr += ' # text={{Master_Id={}}}\n'.format(i + 1)
                ofh.write(ostr)

            print("Created: {}".format(outfile))

        if outline['status'] == 'qa':
            dump_qa(ensemble, outdir, i + 1, outline,
                    creator=creator,
                    revision=revision,
                    color=qa_color)

    # Write out the "cmst3" file, as a FITS file.
    #
    filename = 'master_hulls.{}.v{:03d}.fits'.format(ensemble,
                                                     revision)
    outfile = os.path.join(outdir, filename)
    write_hulls(ensemble, outfile, master_hulls, hull_areas,
                outlines,
                hull_store,
                stacks=stacks,
                compzero=compzero,
                revision=revision,
                creator=creator)


if __name__ == "__main__":

    import argparse
    import sys

    parser = argparse.ArgumentParser(description=help_str,
                                     prog=sys.argv[0])

    parser.add_argument("ensemblefile",
                        help="The ensemble to stack mapping")
    parser.add_argument("ensemble", type=str,
                        help="The ensemble to process")
    parser.add_argument("outdir", type=str,
                        help="The output directory (must not exist)")

    parser.add_argument("--mrgsrc3dir",
                        default="/data/L3/chs_master_match/input/mrgsrc3",
                        help="The mrgsrc3 directory: default %(default)s")
    parser.add_argument("--compzero", type=int,
                        default=9000,
                        help="The COMPONENT value for stack hull cpt=0: default %(default)s")

    args = parser.parse_args(sys.argv[1:])

    process_ensemble(args.ensemblefile, args.ensemble, args.outdir,
                     mrgsrc3dir=args.mrgsrc3dir,
                     compzero=args.compzero,
                     creator=sys.argv[0])
