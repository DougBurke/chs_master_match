"""
Given a list of hull components, create a "merged" hull.

The algorithm is basically

  - create a 1/0 image for each hull (1 for inside the hull)
  - add them up (i.e. they are all on the same grid)
  - smooth them slightly
  - create a contour at a fraction of the number of stacks
    containing hulls (not the number of hulls, as can have multiple
    hulls from the same stack in a master hull)

"""

import os
import subprocess
import tempfile

import six

import numpy as np

import cxcdm
import paramio
import pycrates

import chs_utils as utils


# --- START: Kenny's convex-hull code ---

# Lifted from : http://code.activestate.com/recipes/66527-finding-the-convex-hull-of-a-set-of-2d-points/
# Replaced with numpy's determinant which is more numerically stable

class ConvexHull():

    def __init__(self, x, y):
        self.hull = self.convexHull(zip(x, y))
        self.cx = [xy[0] for xy in self.hull]
        self.cy = [xy[1] for xy in self.hull]

    @staticmethod
    def _myDet(p, q, r):
        """Calc. determinant of a special matrix with three 2D points.

        The sign, "-" or "+", determines the side, right or left,
        respectivly, on which the point r lies, when measured against
        a directed vector from p to q.
        """

        a = np.array([[1.0, p[0], p[1]],
                      [1.0, q[0], q[1]],
                      [1.0, r[0], r[1]]])
        det = np.linalg.det(a)
        return det

    def _isRightTurn(self, z):
        "Do the vectors pq:qr form a right turn, or not?"
        (p, q, r) = z
        assert p != q and q != r and p != r

        if self._myDet(p, q, r) < 0:
            return 1
        else:
            return 0

    def convexHull(self, P):
        "Calculate the convex hull of a set of points."

        # Get a local list copy of the points and sort them lexically.
        # points = map(None, P)
        points = list(P)
        points.sort()

        # Build upper half of the hull.
        upper = [points[0], points[1]]
        for p in points[2:]:
            upper.append(p)
            while len(upper) > 2 and not self._isRightTurn(upper[-3:]):
                del upper[-2]

        # Build lower half of the hull.
        points.reverse()
        lower = [points[0], points[1]]
        for p in points[2:]:
            lower.append(p)
            while len(lower) > 2 and not self._isRightTurn(lower[-3:]):
                del lower[-2]

        # Remove duplicates.
        del lower[0]
        del lower[-1]

        # Concatenate both halfs and return.
        return tuple(upper + lower)


# --- END: Kenny's convex-hull code ---

# --- START: Kenny's co-linear code ---

def remove_colinear(xy, threshold):
    """
    This isn't a great algorithm but does a pretty good
    first-pass kind of job

    Computes the area of triangle for 3 points. If area is too big
    then we keep the middle point, otherwise we move onto the next
    one.
    """
    retval = [xy[0]]  # always keep the 1st point
    for ii in six.moves.xrange(2, len(xy)):  # Polygon must have 3 points!

        x0 = retval[-1][0]  # The last saved point
        y0 = retval[-1][1]

        x2 = xy[ii][0]  # The current end-point
        y2 = xy[ii][1]

        x1 = xy[ii - 1][0]  # The middle point
        y1 = xy[ii - 1][1]

        area = np.abs(0.5 * ((x1 - x0) * (y2 - y0) -
                             (x2 - x0) * (y1 - y0)))
        if area > threshold:
            # If area is too big, then keep last point
            retval.append(xy[ii - 1])

    return retval


# --- END: Kenny's co-linear code ---


# Note: the boundary handling could have been done as an object
#       but I didn't feel it was worth it.
#
def find_poly_bounds(poly):
    """Return the bounds of the polygon.

    Parameters
    ----------
    poly : numpy array
        The polygon data, with a shape of (2, npts). All values
        are expected to be finite.

    Returns
    -------
    bounds : dict
        The bounds of the polygon with the keys
        xmin, xmax, ymin, ymax.
    """

    x = poly[0, :]
    y = poly[1, :]
    return {'xmin': x.min(),
            'xmax': x.max(),
            'ymin': y.min(),
            'ymax': y.max()}


def empty_bounds():
    """Return an empty bounds dictionary."""

    return {'xmin': None, 'xmax': None, 'ymin': None, 'ymax': None}


def validate_bounds(bounds):
    """Error if any element is not set or invalid."""

    for n in ['xmin', 'ymin', 'xmax', 'ymax']:
        if bounds[n] is None:
            raise ValueError("{} key is not set".format(n))

    for n in ['x', 'y']:
        lo = bounds[n + 'min']
        hi = bounds[n + 'max']
        # allow equality
        if lo > hi:
            raise ValueError("{} axis is not ordered".format(n))


def normalize_bounds(bounds, pixscale=1):
    """Make bounds fall on 'nice' pixel boundaries.

    Parameters
    ----------
    bounds : dict
        Must have the keywords 'xmin', 'xmax', 'ymin', and
        'ymax'. These values may be updated, and they can not
        be None.
    pixscale : int, optional
        The binning factor used for the output image, in units
        of the native pixel scale. It is assumed that pixscale is
        1 or greater.


    Notes
    -----
    The idea is to return limits which are on a "native"
    scale: that is, pixel coordinates are set up so that
    when native pixel binning is used the edges are n-0.5 to
    n+0.5 for pixel n (n is an integer).

    Ideally it would be tweaked so that - when binning is used -
    the pixels are still aligned with the 0.5 edge.

    A pixel value of p + dp where dp = [0.5, 1.5) and p is an
    integer will be mapped to p + 0.5 (minimum) and p + 1.5
    for maximum.

    The maximum ranges are tweaked so that they cover the full
    pixel (only relevant when pixscale > 1).
    """

    assert pixscale >= 1
    assert int(pixscale) == pixscale

    for n in ['xmin', 'ymin', 'xmax', 'ymax']:
        if bounds[n] is None:
            raise ValueError("{} key is not set".format(n))

    # This should work if the bounds are negative, but has seen
    # limited testing.
    #
    for n in ['xmin', 'ymin']:
        pval = np.floor(bounds[n] + 0.5)
        bounds[n] = pval - 0.5

    for n in ['xmax', 'ymax']:
        pval = np.floor(bounds[n] + 0.5)
        bounds[n] = pval + 0.5

    # Now check that the limits work with the requested pixel scale
    #
    nx = bounds['xmax'] - bounds['xmin']
    ny = bounds['ymax'] - bounds['ymin']

    assert int(nx) == nx
    assert int(ny) == ny

    nx = int(nx)
    ny = int(ny)
    pixscale = int(pixscale)

    dx = nx % pixscale
    if dx > 0:
        bounds['xmax'] += (pixscale - dx)
        nx += (pixscale - dx)

    dy = ny % pixscale
    if dy > 0:
        bounds['ymax'] += (pixscale - dy)
        ny += (pixscale - dy)

    assert nx % pixscale == 0
    assert ny % pixscale == 0

    nx = nx // pixscale
    ny = ny // pixscale

    bounds['nx'] = nx
    bounds['ny'] = ny
    bounds['pixsize'] = pixscale


def expand_bounds(bounds, delta):
    """Expand the boundary by delta image pixels.

    Both the upper and lower bounds are changed, so the
    width and height of the rectangle change by 2 * delta.
    This should be called before normalize_bounds is called,
    particularly if pixscale is > 1.

    Parameters
    ----------
    bounds : dict
        Must have the keywords 'xmin', 'xmax', 'ymin', and
        'ymax'. These values may be updated, and they can not
        be None.
    delta : int
        The number of pixels (image) to expand the boundary by;
        if <= 0 nothing is done.
    """

    if delta <= 0:
        return

    for n in ['xmin', 'ymin', 'xmax', 'ymax']:
        if bounds[n] is None:
            raise ValueError("{} key is not set".format(n))

    for n in ['xmin', 'ymin']:
        bounds[n] -= delta

    for n in ['xmax', 'ymax']:
        bounds[n] += delta


def update_bounds(bounds, newbounds):
    """Update the bounds if necessary.

    Parameters
    ----------
    bounds, newbounds : dict
        Must have the keywords 'xmin', 'xmax', 'ymin', and
        'ymax'. Note that bounds may be updated by this routine.
        It is required that newbounds has concrete values, but
        bounds may have None values.
    """

    for n in ['xmin', 'ymin']:
        assert newbounds[n] is not None
        if bounds[n] is None or bounds[n] > newbounds[n]:
            bounds[n] = newbounds[n]

    for n in ['xmax', 'ymax']:
        assert newbounds[n] is not None
        if bounds[n] is None or bounds[n] < newbounds[n]:
            bounds[n] = newbounds[n]


def make_ones_image(outfile, bounds, transform):
    """Create an image of all 1's.

    Parameters
    ----------
    outfile : str
        The name of the image to create. It will be clobbered if
        it already exists.
    bounds : dict
        The bounds of the image with the keys
        xmin, xmax, ymin, ymax, nx, ny, and pixsize. It is assumed to
        have been passed through normalize_bounds() and that each
        dimension has at least 1 pixel.
    transform
        The SKY to celestial transformation for the image. Assumed to
        be a TANGENT transform.

    Notes
    -----
    The image has a datatype of int16 / Int2 as this seems to be
    the easiest to create. uint8 should be all that's needed but
    memory should not be a problem here.
    """

    ctype = transform.get_parameter_value('CTYPE')
    assert ctype[0] == 'RA---TAN'
    assert ctype[1] == 'DEC--TAN'

    nx = bounds['xmax'] - bounds['xmin']
    ny = bounds['ymax'] - bounds['ymin']
    assert int(nx) == nx
    assert int(ny) == ny
    assert nx > 0
    assert ny > 0

    nx = bounds['nx']
    ny = bounds['ny']
    assert int(nx) == nx
    assert int(ny) == ny
    assert nx > 0
    assert ny > 0

    pixsize = bounds['pixsize']
    assert pixsize >= 1
    assert int(pixsize) == pixsize

    isize = np.asarray([ny, nx], dtype=np.int32)

    if os.path.exists(outfile):
        os.remove(outfile)

    # Seems like 16-bit integer is easiest to create
    itype = np.int16

    # It is easier to do this via cxcdm than crates at this time
    ds = cxcdm.dmDatasetCreate(outfile)
    bl = cxcdm.dmDatasetCreateImage(ds, 'COVERAGE', itype, isize)
    dd = cxcdm.dmImageGetDataDescriptor(bl)

    # note that the bounds values give the edges of the pixel,
    # so it makes sense to use crpix = 0.5 for the low value.
    #
    sky = cxcdm.dmArrayCreateAxisGroup(dd, 'SKY', np.float64,
                                       'pixel', ['X', 'Y'])
    crpix = np.asarray([0.5, 0.5])
    crval = np.asarray([bounds['xmin'], bounds['ymin']])
    cdelt = np.asarray([pixsize, pixsize])
    cxcdm.dmCoordSetTransform(sky, crpix, crval, cdelt)

    # NOTE: even if the image is binned, the CDELT value should
    # not need changing (as it relates to SKY pixels).
    #
    crpix = transform.get_parameter_value('CRPIX')
    crval = transform.get_parameter_value('CRVAL')
    cdelt = transform.get_parameter_value('CDELT')

    # Note: do nothing with the descriptor that is returned
    #       and we do not set the "extra" parameters argument
    #       to avoid creating LONGPOLE/LATPOLE keywords
    cxcdm.dmCoordCreate(sky, 'EQPOS', 'degree',
                        ['RA', 'DEC'], 'TAN',
                        crpix, crval, cdelt)

    # set all pixels to 1
    ivals = np.ones(isize, dtype=itype)
    cxcdm.dmImageSetData(dd, ivals)

    cxcdm.dmDatasetClose(ds)


def get_region_filters(infile, hulls, ndigit=2):
    """Create list of spatially-filtered versions of infile.

    Parameters
    ----------
    infile : str
        The input image file name. The file is assumed to have a
        SKY vector coordinate in the correct system.
    hulls : list of polygons
        The hulls to use to filter the input; each element of
        the list is assumed to be a 2D NumPy array, wish shape
        (2,npts), ans in the SKY coordinates of the input image.
    ndigit : int, optional
        The number of decimal places to use when converting the
        polygons to a filter. If None then let Python decide.
        There is no check to see if the polygon remains simple
        after this change to the coordinate values.

    Returns
    -------
    fnames : list of str
        The list of file names, including spatial filters representing
        the hulls. There is no check on the length of these names.
        The files include an '[opt full]' specifier so that the image
        sizes are unchanged.
    """

    fnames = []
    for poly in hulls:
        rstr = utils.make_region_string(poly, ndp=ndigit)
        fnames.append("{}[sky={}][opt full]".format(infile, rstr))

    return fnames


def combine_filtered_images(infile, outfile, hulls,
                            ndigit=2, maxlen=4000,
                            tmpdir=None):
    """Filter infile by the hulls and add up.

    Parameters
    ----------
    infile : str
        The input image (created by make_ones_image). It is assumed
        to exceed the bounds of all the hulls, but this assumption
        is not needed here.
    outfile : str
        The image to create (will be clobbered).
    hulls : list of polygons
        The hulls to use to filter the input; each element of
        the list is assumed to be a 2D NumPy array, with shape
        (2, npts), and in the SKY coordinates of the input image.
    ndigit : int, optional
        The number of decimal places to use when converting the
        polygons to a filter. If None then let Python decide.
        There is no check to see if the polygon remains simple
        after this change to the coordinate values.
    maxlen : int, optional
        Warn if any file name (that is infile + spatial filter
        expression) exceeds this length (if not None). The default
        value is somewhat made up, and given that a stack file
        is used for the filename inputs to dmimgcalc may not be
        needed.
    tmpdir : None or str
        The directory to use for temporary files. If None, the default
        is used (the Python tempfile isn't very explicit about
        what this default is, but it tends to be /tmp on Linux).

    """

    nhulls = len(hulls)
    assert nhulls > 1

    fnames = get_region_filters(infile, hulls, ndigit=ndigit)

    # quick check on the length of these file names, just in case
    if maxlen is not None:
        assert maxlen > 0
        for i, fname in enumerate(fnames):
            if len(fname) > maxlen:
                print("WARNING: maxlen exceeded for hull " +
                      "#{}".format(i + 1))

    # Run dmimgcalc: use a stack for the input files
    #
    stackfile = tempfile.NamedTemporaryFile(prefix='hulltmp',
                                            suffix='.stk',
                                            dir=tmpdir)
    with open(stackfile.name, 'w') as fh:
        for fname in fnames:
            fh.write(fname + "\n")

    # add up all the images
    imgsum = "+".join(['img{}'.format(i + 1) for i in range(nhulls)])

    cmd = ['dmimgcalc',
           'infile=@' + stackfile.name,
           'infile2=',
           'outfile=' + outfile,
           'operation=imgout=' + imgsum,
           'mode=h',
           'clobber=yes']

    paramio.punlearn(cmd[0])
    subprocess.check_call(cmd)


def smooth_image(infile, outfile, sigma=3, nsigma=5):
    """Gaussian smooth the image.

    Note that the convolution is done using a sliding-cell rather
    than FFT because: it is easier to ensure no edge effects
    since can set edges=constant when method=slide, and this is
    important for running dmcontour on the output; I found that
    numerical issues from the FFT appeared to lead to more
    'problem shapes' created by dmcontour; image sizes do not
    get padded.

    Parameters
    ----------
    infile : str
        The input image.
    outfile : str
        The image to create (will be clobbered).
    sigma : float, optional
        The gaussian width, in pixels.
    nsigma : float, optional
        The size of the box, in units of sigma.
    """

    kernel = "lib:gaus(2,{0},1,{1},{1})".format(nsigma, sigma)

    cmd = ['aconvolve',
           'infile=' + infile,
           'outfile=' + outfile,
           'kernelspec=' + kernel,
           'method=slide',
           'edges=constant',
           'const=0',
           'clobber=yes']

    paramio.punlearn(cmd[0])
    subprocess.check_call(cmd)


# When using SKY coordinate systems and the ACIS scale,
# I have found that an "area threshold" of 1 gives good results
# for the remove_colinear step. I have not tried to fine tune
# this value.
#
_linear_threshold = 1.0


def make_convex_hull(poly,
                     linear=_linear_threshold):
    """Create a convex hull.

    This should really be SKY coordinates, rather than celestial,
    since there is an area calculation (as a validation step) and
    the threshold for removing colinear points assumes that the
    coordinates are in ACIS pixels (it should be okay to use
    HRC pixels).

    Parameters
    ----------
    poly
        The polygon to turn into a convex hull, with shape (2, npts),
        and all values are finite.
    linear : float, optional
        This is used in the co-linear check, and represents an area.

    Returns
    -------
    hull
        A convex hull, with shape (2, npts), although note that this
        value of npts need not match the input values.

    Notes
    -----
    The coordinates are passed through a simple filter first, to
    remove "co-linear" points - as it makes the convex-hull routine
    more reliable (at least for the data I have used). Note that the
    dmcontour output, which this is expected to be used on, can have
    many points along what is essentially a linear edge, so this step
    is particularly beneficial.
    """

    x = poly[0]
    y = poly[1]
    xy = list(zip(x, y))

    nxy = remove_colinear(xy, linear)
    nx, ny = list(zip(*nxy))

    nx = np.asarray(nx)
    ny = np.asarray(ny)

    hull = ConvexHull(nx, ny)
    hx = np.asarray(hull.cx)
    hy = np.asarray(hull.cy)

    hpoly = np.vstack((hx, hy))
    return utils.validate_polygon(hpoly, report=False)


def contour_image(infile, level, tmpdir=None):
    """Contour an image to create a region.

    If a single polygon is returned then it is guaranteed to be a
    convex hull. This routine does not return the original data
    (i.e. the dmcontour output before converting to a convex hull;
    this information may be useful for QA). It is not yet clear what
    is best to return when multiple polygons need to be returned,
    including the case of excluded polygons.

    Parameters
    ----------
    infile : str
        The image to contour.
    level : float
        The contour level
    tmpdir : None or str
        The directory to use for temporary files. If None, the default
        is used (the Python tempfile isn't very explicit about
        what this default is, but it tends to be /tmp on Linux).

    Returns
    -------
    out : dict
        The status keyword is one of: 'okay', 'qa', or 'failed'.
        If status is 'okay' then the convex hull is returned in the
        'hull' key, and is in celestial coordinates.
        The reason keyword is set if status is not 'okay' and provides
        more information; other keywords may be set in certain
        situations to return more data: the 'hulls' keyword contains
        a list of (shape, coordinate) pairs and can be the empty
        list: the coordinates are in the celestial system, and
        when status='qa' the values are convex hulls for included
        polygons. The logic here needs tweaking for the pipeline,
        so I have not spent much time coming up with anything
        sensible.

        TODO: this documention has not been updated to match the
        current behavior

    """

    outfile = tempfile.NamedTemporaryFile(prefix='hulltmp',
                                          suffix='.contour',
                                          dir=tmpdir)

    cmd = ['dmcontour',
           'infile=' + infile,
           'outfile=' + outfile.name,
           'levels={}'.format(level),
           'clobber=yes']

    # It is possible that dmcontour can fail, so we need to handle
    # this case.
    #
    paramio.punlearn(cmd[0])
    try:
        subprocess.check_call(cmd)
    except subprocess.CalledProcessError:
        # It is not clear whether is is worth logging the error
        #
        return {'status': 'failed',
                'reason': 'dmcontour failed'}

    # NOTE: there is no filter on the SHAPE field since there are
    # problems in CIAO 4.9 (and earlier) with string filtering,
    # so it is safest to filter in Python
    #
    cr = pycrates.read_file(outfile.name)
    nrows = cr.get_nrows()
    if nrows == 0:
        return {'status': 'failed',
                'reason': 'dmcontour returned an empty file'}

    # I am assuming that the only shape values created by
    # dmcontour are Polygon, !Polygon, or empty (which happens
    # when there has been an error but dmcontour in CIAO 4.8 hasn't
    # handled it).
    #
    shapes = cr.SHAPE.values
    idx = np.where(shapes != '')
    shapes = shapes[idx]

    # Convert from polygon to convex hull in the SKY coordinate
    # system, as this is probably safer than doing it to celestial
    # coordinates (may tests show little difference, but the SKY
    # system should be better for edge cases such as if the source
    # crosses ra=0/360 degrees or covers a pole).
    #
    # Assume the names are POS for SKY and EQPOS for SKY to
    # celestial (needed since the return value is in celestial
    # coordinates).
    #
    tr = cr.get_transform('EQPOS')

    def sky_to_cel(poly):
        # convert from (2,npts) to (1,2,npts) to
        # make apply happy, then strip off the extra
        # dimension
        invals = poly[np.newaxis]
        outvals = tr.apply(invals)[0]
        return outvals

    sky = []
    for poly in cr.get_column('POS').values[idx]:
        cleaned = utils.validate_polygon(poly, report=False)
        sky.append(cleaned)

    npoly = len(sky)
    if npoly != nrows:
        # return the data that was created, if any
        #
        # The coordinates *could* be turned into convex hulls
        # here, but for now leave as is.
        #
        celestial = []
        for poly in sky:
            celestial.append(sky_to_cel(poly))

        return {'status': 'failed',
                'reason': 'dmcontour returned 1 or more blank shapes',
                'hulls': list(zip(shapes, celestial))}

    # If there is one row of a Polygon then convert to a convex
    # hull and return.
    if npoly == 1:
        sky_hull = make_convex_hull(sky[0])
        cel_hull = sky_to_cel(sky_hull)
        return {'status': 'okay',
                'hull_cel': cel_hull,
                'hull_sky': sky_hull}

    # Special case 1 included and n excluded polygons.
    # This has been shoe-horned into the code and could probably
    # be a lot neater.
    #
    scount = {}
    for shape in shapes:
        try:
            scount[shape] += 1
        except KeyError:
            scount[shape] = 1

    if len(scount) > 2 or 'Polygon' not in scount:
        return {'status': 'error',
                'reason': 'unexpected shapes',
                'hulls': list(zip(shapes, sky))}

    if scount['Polygon'] == 1:
        polys = [poly for shape, poly in zip(shapes, sky)
                 if shape == 'Polygon']
        assert len(polys) == 1
        sky_hull = make_convex_hull(polys[0])
        cel_hull = sky_to_cel(sky_hull)
        return {'status': 'okay',
                'hull_cel': cel_hull,
                'hull_sky': sky_hull
                }

    # If there are multiple hulls, note as a qa case, and return
    # the values for the INCLUDED polygons only.
    #
    hulls_cel = []
    hulls_sky = []
    for (shape, poly) in zip(shapes, sky):
        if shape == 'Polygon':
            sky_hull = make_convex_hull(poly)
            cel_hull = sky_to_cel(sky_hull)
        elif shape == '!Polygon':
            # for now just skip
            continue
        else:
            raise RuntimeError("should not be possible")

        hulls_cel.append(cel_hull)
        hulls_sky.append(sky_hull)

    return {'status': 'qa',
            'reason': 'multiple hulls found',
            'hulls_cel': hulls_cel,
            'hulls_sky': hulls_sky}


def read_hull(stkid, cpt, indir):
    """Read in the convex hull data for the stack.

    Parameters
    ----------
    stkid : str
        The stack id.
    cpt : int
        The component number
    indir : str
        The location of the mrgsrc3 file.

    """

    mrg3file = utils.find_mrgsrc3(stkid, indir)
    infile = "{}[MEXTSRC][STATUS=0,COMPONENT={}]".format(mrg3file,
                                                         cpt)
    cr = pycrates.read_file(infile)
    if cr.get_nrows() != 1:
        raise IOError("No valid hull for {} {}".format(stkid, cpt))

    pos = cr.POS.values.copy()[0]
    eqsrc = cr.get_column('eqsrc').values.copy()[0]
    tr = cr.get_transform('EQSRC').copy()

    # remove any non-finite values (assume same in POS and EQSRC,
    # and X/Y).
    #
    idx = np.isfinite(pos[0])
    pos = pos[:, idx]
    eqsrc = eqsrc[:, idx]

    return {'cpt': cpt,
            'tr': tr,
            'pos': pos,
            'eqsrc': eqsrc,
            'stkid': stkid,
            'infile': infile}


def merge_polygons(hulls,
                   acceptfrac=0.2,
                   maxcount="cohorts",
                   sigma=3,
                   nsigma=5,
                   tmpdir=None
                   ):
    """Merge overlapping polygons.

    Parameters
    ----------
    hulls : list of dict
        The hulls to merge; these are the output of read_hull
    acceptfrac : float, optional
        The fraction at which to draw the merged polygon
        (when multiple hulls are present). The value used
        is acceptfrac * n, where n is determined by maxcount.
    maxcount : {'cohorts', 'hulls'}
        Should n be the number of different cohorts in the
        list of overlapping hulls (maxcount='cohorts', or the
        number of convex hulls (maxcount='hulls'). This only
        makes a difference when there are multiple hulls from
        a single cohort in the list.
    nsigma, sigma : float
        The number of sigma (the box size) and the sigma, in pixels,
        of the gaussian used to smooth the image. If either is None
        then no smoothing is done. Note that the smoothing scale,
        sigma, is given in ACIS pixels, so is multiplied by 3.8
        before being applied to HRC data (so that the physical
        scale being smoothed is similar).
    tmpdir : None or str
        The directory to use for temporary files. If None, the default
        is used (the Python tempfile isn't very explicit about
        what this default is, but it tends to be /tmp on Linux).

    Returns
    -------
    res, stkid : (dict, str)
        res is a dictionary listing the results: see contour_image for
        a brief discussion of the keys although it also has the
        'inputs' key, which contains the polygons from
        the input sources (in celestial coordinates). It is a list
        of polygons (so NumPy arrays with shape (2, npts), where
        the npts value is not necessarily the same for every
        polygon). stkid is a string, and lists the STACKID of the
        stack used for the calculation; it is intended that this is
        a repeatable value (i.e. re-running will get the same stack,
        but uses a rather simple scheme).

    """

    assert len(hulls) > 0

    if acceptfrac <= 0 or acceptfrac > 1:
        raise ValueError("acceptfreac should be (0,1]")

    smoothing = not (sigma is None or nsigma is None)
    if smoothing:
        if sigma <= 0:
            raise ValueError("sigma > 0, not {}".format(sigma))
        if nsigma <= 0:
            raise ValueError("nsigma > 0, not {}".format(nsigma))

    # What is the value that we apply acceptfrac to?
    #
    if maxcount == "cohorts":
        cohorts = set([h['stkid'] for h in hulls])
        maxnum = len(cohorts)

    elif maxcount == "hulls":
        maxnum = len(hulls)

    else:
        raise ValueError("invalid maxcount='{}'".format(maxcount))

    # Need to pick a base coordinate system to do the calculations
    # in. There is no real reason to pick any one in particular
    # (the closest tangent point to this source may be nice but
    # not worth the effort). I am going to pick the first ACIS
    # cohort I find (falling back to HRC if necessary).
    #
    # This is not the most-elegant pieces of code, but I doubt
    # it is a major time sink.
    #
    # The hull list is sorted so that the choice of base stack is
    # repeatable.
    #
    tr = None
    tr_stackid = None
    is_acis = False

    for h in sorted(hulls, key=lambda h: h['stkid']):
        if h['stkid'].startswith('acis'):
            tr = h['tr']
            tr_stackid = h['stkid']
            is_acis = True
            break
        elif tr is None:
            tr = h['tr']
            tr_stackid = h['stkid']

    # Scale the smoothing scale so that it matches the HRC pixel
    # size, if necessary. That is, we want to smooth by a kernel
    # with approximately the same size in arcsec whether it is
    # ACIS or HRC data.
    if sigma is not None and not is_acis:
        sigma *= 3.8

    # Convert polygons to SKY coordinates and calculate the bounds
    # of each polygon in sky coordinates. The SKY coordinate system
    # is from the "base cohort" chosen above.
    #
    polys = []
    bounds = empty_bounds()
    for hull in hulls:

        # coords is (2, npts) and it is easier to convert if
        # (1, 2, npts).
        #
        coords = hull['eqsrc'][np.newaxis]
        sky = tr.invert(coords)[0]
        polys.append(sky)

        pbounds = find_poly_bounds(sky)
        update_bounds(bounds, pbounds)

    # Calculate the bounds for the image. We add on a boundary
    # to make sure that dmcontour does not get confused (we do
    # not want any contour to hit the edge of the image, so
    # we want a few pixels of empty space at each edge).
    # Kenny has suggested a fix to dmcontour to address this,
    # but I am assuming it is not available.
    #
    # The image coordinates do not need to fall onto "nice"
    # pixel boundaries, but it's easier for my testing and
    # verification.
    #
    validate_bounds(bounds)

    boundary = 2
    if smoothing:
        # I think that nsigma/2 could be used here, but it is
        # not entirely clear, and the extra space doesn't have
        # a huge impact on runtime/memory use.
        #
        boundary += np.ceil(nsigma * sigma)

    expand_bounds(bounds, boundary)
    normalize_bounds(bounds, pixscale=1)

    # Create a temporary image which is all 1's that covers
    # this boundary range.
    #
    basefile = tempfile.NamedTemporaryFile(prefix='hulltmp',
                                           suffix='.img',
                                           dir=tmpdir)
    make_ones_image(basefile.name, bounds, tr)

    # Apply each polygon to this file and sum up.
    #
    filtfile = tempfile.NamedTemporaryFile(prefix='hulltmp',
                                           suffix='.fimg',
                                           dir=tmpdir)
    combine_filtered_images(basefile.name, filtfile.name,
                            polys, tmpdir=tmpdir)

    # Smooth the image
    #
    if smoothing:
        smoothfile = tempfile.NamedTemporaryFile(prefix='hulltmp',
                                                 suffix='.simg',
                                                 dir=tmpdir)
        smooth_image(filtfile.name, smoothfile.name,
                     sigma=sigma, nsigma=nsigma)
        filtfile = smoothfile

    # Create outline.
    #
    level = acceptfrac * maxnum
    out = contour_image(filtfile.name, level,
                        tmpdir=tmpdir)

    return out, tr_stackid


def make_merged_hull(hullcpts, mrgsrc3dir,
                     acceptfrac=0.2,
                     maxcount="cohorts",
                     sigma=3, nsigma=5,
                     tmpdir=None):
    """Create the merged hull.

    Parameters
    ----------
    hullcpts : sequence of (stackid, component)
        The stack-level hulls to merge
    mrgsrc3dir : str
        The location of the mrgsrc3 files.
    acceptfrac : float, optional
        The fraction at which to draw the merged polygon
        (when multiple hulls are present). The value used
        is acceptfrac * n, where n is determined by maxcount.
    maxcount : {'cohorts', 'hulls'}
        Should n be the number of different cohorts in the
        list of overlapping hulls (maxcount='cohorts', or the
        number of convex hulls (maxcount='hulls'). This only
        makes a difference when there are multiple hulls from
        a single cohort in the list.
    nsigma, sigma : float
        The number of sigma (the box size) and the sigma, in pixels,
        of the gaussian used to smooth the image. If either is None
        then no smoothing is done. Note that the smoothing scale,
        sigma, is given in ACIS pixels, so is multiplied by 3.8
        before being applied to HRC data (so that the physical
        scale being smoothed is similar).
    tmpdir : None or str
        The directory to use for temporary files. If None, the default
        is used (the Python tempfile isn't very explicit about
        what this default is, but it tends to be /tmp on Linux).

    Returns
    -------
    outline : dict
        The outline of the merged hull. The keys are:
        status, eqpos, pos, base_stack. If status is not 'okay'
        then eqpos and pos can be None (status is 'error', in
        which case base_stack will also be None)
        or 3D shapes (i.e. npolygons, 2, npts), when status is
        'qa'.

    Notes
    -----
    Information flow on failure (including unexpected or qa) could
    be better.
    """

    if len(hullcpts) == 0:
        raise ValueError("hullcpts is empty")

    hulls = []
    for stkid, cpt in hullcpts:
        hulls.append(read_hull(stkid, cpt, mrgsrc3dir))

    # If there is only one hull then just promote the stack-level
    # data.
    #
    if len(hulls) == 1:
        h0 = hulls[0]
        outline = {'status': 'okay',
                   'eqpos': h0['eqsrc'],
                   'pos': h0['pos'],
                   'base_stack': h0['stkid']}
        return outline

    out, tr_stkid = merge_polygons(hulls, acceptfrac=0.2,
                                   maxcount="cohorts",
                                   sigma=3, nsigma=5,
                                   tmpdir=None)

    # ARGH: the keys are different to read_hull, which makes this
    #       needlessly different, but not worth changing now.
    #
    if out['status'] == 'okay':
        outline = {'status': 'okay',
                   'eqpos': out['hull_cel'],
                   'pos': out['hull_sky'],
                   'base_stack': tr_stkid}
        return outline

    print("WARNING: when calculating merged hull got")
    print("status= {}".format(out['status']))
    print("reason= {}".format(out['reason']))

    # Possible failures are
    #    status=failed   - may not have any data
    #    status=error    - returned unexpected data
    #    status=qa       - multiple hulls returned
    #
    if out['status'] == 'qa':
        outline = {'status': 'qa',
                   'eqpos': out['hulls_cel'],
                   'pos': out['hulls_sky'],
                   'base_stack': tr_stkid}
        return outline

    # There could be some data to return here, but for now
    # error out everything.
    return {'status': 'error',
            'eqpos': None,
            'pos': None,
            'base_stack': None}


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
