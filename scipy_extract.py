# this is a dropin replacement for the following scipy funcions:
# zoom, median_filter and gaussian_filter
#
# instead of importing from scipy as:
#   from scipy.ndimage import zoom 
#   from scipy.ndimage import median_filter 
#   from scipy.ndimage import gaussian_filter 
#   from scipy.ndimage import label
#   from scipy.ndimage.interpolation import rotate
# you can import from here as:   
#   from scipy_extract import zoom
#   from scipy_extract import median_filter
#   from scipy_extract import gaussian_filter
#   from scipy_extract import label
#   from scipy_extract import rotate
#
# be aware that you will also need the files
#    scipy\ndimage\_nd_image.pyd
# and
#    scipy\ndimage\_ni_label.pyd
# and
#    scipy\optimize\_minpack.pyd
#    scipy\extra-dll\libchkder.6HLXPVTQJEGRZGLI5DFRMNW3SS76BHP6.gfortran-win_amd64.dll
# to be copied to the directory of this script
#

import numpy as np
import sys
import math
import threading
import _nd_image
import _ni_label
import _minpack

# --- extracted from:  scipy/_lib/six.py ---

if sys.version_info[0] == 3:
    string_types = str,
else:
    string_types = basestring,

    
# --- extracted from:  scipy/ndimage/_ni_support.py ---

def _normalize_sequence(input, rank, array_type=None):
    is_str = isinstance(input, string_types)
    if hasattr(input, '__iter__') and not is_str:
        normalized = list(input)
        if len(normalized) != rank:
            err = "sequence argument must have length equal to input rank"
            raise RuntimeError(err)
    else:
        normalized = [input] * rank
    return normalized

def _get_output(output, input, shape=None):
    if shape is None:
        shape = input.shape
    if output is None:
        output = np.zeros(shape, dtype=input.dtype.name)
    elif type(output) in [type(type), type(np.zeros((4,)).dtype)]:
        output = np.zeros(shape, dtype=output)
    elif type(output) in string_types:
        output = np.typeDict[output]
        output = np.zeros(shape, dtype=output)
    elif output.shape != shape:
        raise RuntimeError("output shape not correct")
    return output    

def _extend_mode_to_code(mode):
    if mode == 'nearest':
        return 0
    elif mode == 'wrap':
        return 1
    elif mode == 'reflect':
        return 2
    elif mode == 'mirror':
        return 3
    elif mode == 'constant':
        return 4
    else:
        raise RuntimeError('boundary mode not supported')    

def _check_axis(axis, rank):
    if axis < 0:
        axis += rank
    if axis < 0 or axis >= rank:
        raise ValueError('invalid axis')
    return axis


# --- extracted from:  scipy/ndimage/filters.py ---
    
def _rank_filter(input, rank, size=None, footprint=None, output=None,
     mode="reflect", cval=0.0, origin=0, operation='rank'):
    input = np.asarray(input)
    if np.iscomplexobj(input):
        raise TypeError('Complex type not supported')
    origins = _normalize_sequence(origin, input.ndim)
    if footprint is None:
        if size is None:
            raise RuntimeError("no footprint or filter size provided")
        sizes = _normalize_sequence(size, input.ndim)
        footprint = np.ones(sizes, dtype=bool)
    else:
        footprint = np.asarray(footprint, dtype=bool)
    fshape = [ii for ii in footprint.shape if ii > 0]
    if len(fshape) != input.ndim:
        raise RuntimeError('filter footprint array has incorrect shape.')
    for origin, lenf in zip(origins, fshape):
        if (lenf // 2 + origin < 0) or (lenf // 2 + origin >= lenf):
            raise ValueError('invalid origin')
    if not footprint.flags.contiguous:
        footprint = footprint.copy()
    filter_size = np.where(footprint, 1, 0).sum()
    if operation == 'median':
        rank = filter_size // 2
    elif operation == 'percentile':
        percentile = rank
        if percentile < 0.0:
            percentile += 100.0
        if percentile < 0 or percentile > 100:
            raise RuntimeError('invalid percentile')
        if percentile == 100.0:
            rank = filter_size - 1
        else:
            rank = int(float(filter_size) * percentile / 100.0)
    if rank < 0:
        rank += filter_size
    if rank < 0 or rank >= filter_size:
        raise RuntimeError('rank not within filter footprint size')
    if rank == 0:
        return minimum_filter(input, None, footprint, output, mode, cval,
                              origins)
    elif rank == filter_size - 1:
        return maximum_filter(input, None, footprint, output, mode, cval,
                              origins)
    else:
        output = _get_output(output, input)
        mode = _extend_mode_to_code(mode)
        _nd_image.rank_filter(input, rank, footprint, output, mode, cval,
                              origins)
        return output

        
def gaussian_filter1d(input, sigma, axis=-1, order=0, output=None,
                      mode="reflect", cval=0.0, truncate=4.0):
    if order not in range(4):
        raise ValueError('Order outside 0..3 not implemented')
    sd = float(sigma)
    # make the radius of the filter equal to truncate standard deviations
    lw = int(truncate * sd + 0.5)
    weights = [0.0] * (2 * lw + 1)
    weights[lw] = 1.0
    sum = 1.0
    sd = sd * sd
    # calculate the kernel:
    for ii in range(1, lw + 1):
        tmp = math.exp(-0.5 * float(ii * ii) / sd)
        weights[lw + ii] = tmp
        weights[lw - ii] = tmp
        sum += 2.0 * tmp
    for ii in range(2 * lw + 1):
        weights[ii] /= sum
    # implement first, second and third order derivatives:
    if order == 1:  # first derivative
        weights[lw] = 0.0
        for ii in range(1, lw + 1):
            x = float(ii)
            tmp = -x / sd * weights[lw + ii]
            weights[lw + ii] = -tmp
            weights[lw - ii] = tmp
    elif order == 2:  # second derivative
        weights[lw] *= -1.0 / sd
        for ii in range(1, lw + 1):
            x = float(ii)
            tmp = (x * x / sd - 1.0) * weights[lw + ii] / sd
            weights[lw + ii] = tmp
            weights[lw - ii] = tmp
    elif order == 3:  # third derivative
        weights[lw] = 0.0
        sd2 = sd * sd
        for ii in range(1, lw + 1):
            x = float(ii)
            tmp = (3.0 - x * x / sd) * x * weights[lw + ii] / sd2
            weights[lw + ii] = -tmp
            weights[lw - ii] = tmp
    return correlate1d(input, weights, axis, output, mode, cval, 0)

def correlate1d(input, weights, axis=-1, output=None, mode="reflect",
                cval=0.0, origin=0):
    input = np.asarray(input)
    if np.iscomplexobj(input):
        raise TypeError('Complex type not supported')
    output = _get_output(output, input)
    weights = np.asarray(weights, dtype=np.float64)
    if weights.ndim != 1 or weights.shape[0] < 1:
        raise RuntimeError('no filter weights given')
    if not weights.flags.contiguous:
        weights = weights.copy()
    axis = _check_axis(axis, input.ndim)
    if (len(weights) // 2 + origin < 0) or (len(weights) // 2 +
                                            origin > len(weights)):
        raise ValueError('invalid origin')
    mode = _extend_mode_to_code(mode)
    _nd_image.correlate1d(input, weights, axis, output, mode, cval,
                          origin)
    return output            
        
def median_filter(input, size=None, footprint=None, output=None,
                  mode="reflect", cval=0.0, origin=0):  
    return _rank_filter(input, 0, size, footprint, output, mode, cval,
                        origin, 'median')  

def gaussian_filter(input, sigma, order=0, output=None,
      mode="reflect", cval=0.0, truncate=4.0):
    input = np.asarray(input)
    output = _get_output(output, input)
    orders = _normalize_sequence(order, input.ndim)
    if not set(orders).issubset(set(range(4))):
        raise ValueError('Order outside 0..4 not implemented')
    sigmas = _normalize_sequence(sigma, input.ndim)
    modes = _normalize_sequence(mode, input.ndim)
    axes = list(range(input.ndim))
    axes = [(axes[ii], sigmas[ii], orders[ii], modes[ii])
            for ii in range(len(axes)) if sigmas[ii] > 1e-15]
    if len(axes) > 0:
        for axis, sigma, order, mode in axes:
            gaussian_filter1d(input, sigma, axis, order, output,
                              mode, cval, truncate)
            input = output
    else:
        output[...] = input[...]
    return output


# --- extracted from:  scipy/ndimage/interpolation.py ---
    
def zoom(input, zoom, output=None, order=3, mode='constant', cval=0.0,
         prefilter=True):
    if order < 0 or order > 5:
        raise RuntimeError('spline order not supported')
    input = np.asarray(input)
    if np.iscomplexobj(input):
        raise TypeError('Complex type not supported')
    if input.ndim < 1:
        raise RuntimeError('input and output rank must be > 0')
    mode = _extend_mode_to_code(mode)
    if prefilter and order > 1:
        filtered = spline_filter(input, order, output=np.float64)
    else:
        filtered = input
    zoom = _normalize_sequence(zoom, input.ndim)
    output_shape = tuple(
            [int(round(ii * jj)) for ii, jj in zip(input.shape, zoom)])

    output_shape_old = tuple(
            [int(ii * jj) for ii, jj in zip(input.shape, zoom)])

    zoom_div = np.array(output_shape, float) - 1
    # Zooming to infinite values is unpredictable, so just choose
    # zoom factor 1 instead
    zoom = np.divide(np.array(input.shape) - 1, zoom_div,
                        out=np.ones_like(input.shape, dtype=np.float64),
                        where=zoom_div != 0)

    output = _get_output(output, input,
                                                   shape=output_shape)
    zoom = np.ascontiguousarray(zoom)
    _nd_image.zoom_shift(filtered, zoom, None, output, order, mode, cval)
    return output

def rotate(input, angle, axes=(1, 0), reshape=True, output=None, order=3,
           mode='constant', cval=0.0, prefilter=True):
    input = np.asarray(input)
    axes = list(axes)
    rank = input.ndim
    if axes[0] < 0:
        axes[0] += rank
    if axes[1] < 0:
        axes[1] += rank
    if axes[0] < 0 or axes[1] < 0 or axes[0] > rank or axes[1] > rank:
        raise RuntimeError('invalid rotation plane specified')
    if axes[0] > axes[1]:
        axes = axes[1], axes[0]
    angle = np.pi / 180 * angle
    m11 = math.cos(angle)
    m12 = math.sin(angle)
    m21 = -math.sin(angle)
    m22 = math.cos(angle)
    matrix = np.array([[m11, m12],
                          [m21, m22]], dtype=np.float64)
    iy = input.shape[axes[0]]
    ix = input.shape[axes[1]]
    if reshape:
        mtrx = np.array([[m11, -m21],
                            [-m12, m22]], dtype=np.float64)
        minc = [0, 0]
        maxc = [0, 0]
        coor = np.dot(mtrx, [0, ix])
        minc, maxc = _minmax(coor, minc, maxc)
        coor = np.dot(mtrx, [iy, 0])
        minc, maxc = _minmax(coor, minc, maxc)
        coor = np.dot(mtrx, [iy, ix])
        minc, maxc = _minmax(coor, minc, maxc)
        oy = int(maxc[0] - minc[0] + 0.5)
        ox = int(maxc[1] - minc[1] + 0.5)
    else:
        oy = input.shape[axes[0]]
        ox = input.shape[axes[1]]
    offset = np.zeros((2,), dtype=np.float64)
    offset[0] = float(oy) / 2.0 - 0.5
    offset[1] = float(ox) / 2.0 - 0.5
    offset = np.dot(matrix, offset)
    tmp = np.zeros((2,), dtype=np.float64)
    tmp[0] = float(iy) / 2.0 - 0.5
    tmp[1] = float(ix) / 2.0 - 0.5
    offset = tmp - offset
    output_shape = list(input.shape)
    output_shape[axes[0]] = oy
    output_shape[axes[1]] = ox
    output_shape = tuple(output_shape)
    output = _get_output(output, input,shape=output_shape)
    if input.ndim <= 2:
        affine_transform(input, matrix, offset, output_shape, output,
                         order, mode, cval, prefilter)
    else:
        coordinates = []
        size = np.product(input.shape, axis=0)
        size //= input.shape[axes[0]]
        size //= input.shape[axes[1]]
        for ii in range(input.ndim):
            if ii not in axes:
                coordinates.append(0)
            else:
                coordinates.append(slice(None, None, None))
        iter_axes = list(range(input.ndim))
        iter_axes.reverse()
        iter_axes.remove(axes[0])
        iter_axes.remove(axes[1])
        os = (output_shape[axes[0]], output_shape[axes[1]])
        for ii in range(size):
            ia = input[tuple(coordinates)]
            oa = output[tuple(coordinates)]
            affine_transform(ia, matrix, offset, os, oa, order, mode,
                             cval, prefilter)
            for jj in iter_axes:
                if coordinates[jj] < input.shape[jj] - 1:
                    coordinates[jj] += 1
                    break
                else:
                    coordinates[jj] = 0
    return output    

def affine_transform(input, matrix, offset=0.0, output_shape=None,
                     output=None, order=3,
                     mode='constant', cval=0.0, prefilter=True):
    if order < 0 or order > 5:
        raise RuntimeError('spline order not supported')
    input = np.asarray(input)
    if np.iscomplexobj(input):
        raise TypeError('Complex type not supported')
    if output_shape is None:
        output_shape = input.shape
    if input.ndim < 1 or len(output_shape) < 1:
        raise RuntimeError('input and output rank must be > 0')
    mode = _extend_mode_to_code(mode)
    if prefilter and order > 1:
        filtered = spline_filter(input, order, output=np.float64)
    else:
        filtered = input
    output = _get_output(output, input,
                                                   shape=output_shape)
    matrix = np.asarray(matrix, dtype=np.float64)
    if matrix.ndim not in [1, 2] or matrix.shape[0] < 1:
        raise RuntimeError('no proper affine matrix provided')
    if (matrix.ndim == 2 and matrix.shape[1] == input.ndim + 1 and
            (matrix.shape[0] in [input.ndim, input.ndim + 1])):
        if matrix.shape[0] == input.ndim + 1:
            exptd = [0] * input.ndim + [1]
            if not np.all(matrix[input.ndim] == exptd):
                msg = ('Expected homogeneous transformation matrix with '
                       'shape %s for image shape %s, but bottom row was '
                       'not equal to %s' % (matrix.shape, input.shape, exptd))
                raise ValueError(msg)
        # assume input is homogeneous coordinate transformation matrix
        offset = matrix[:input.ndim, input.ndim]
        matrix = matrix[:input.ndim, :input.ndim]
    if matrix.shape[0] != input.ndim:
        raise RuntimeError('affine matrix has wrong number of rows')
    if matrix.ndim == 2 and matrix.shape[1] != output.ndim:
        raise RuntimeError('affine matrix has wrong number of columns')
    if not matrix.flags.contiguous:
        matrix = matrix.copy()
    offset = _normalize_sequence(offset, input.ndim)
    offset = np.asarray(offset, dtype=np.float64)
    if offset.ndim != 1 or offset.shape[0] < 1:
        raise RuntimeError('no proper offset provided')
    if not offset.flags.contiguous:
        offset = offset.copy()
    if matrix.ndim == 1:
        _nd_image.zoom_shift(filtered, matrix, offset/matrix, output, order,
                             mode, cval)
    else:
        _nd_image.geometric_transform(filtered, None, None, matrix, offset,
                                      output, order, mode, cval, None, None)
    return output

def spline_filter1d(input, order=3, axis=-1, output=np.float64):
    if order < 0 or order > 5:
        raise RuntimeError('spline order not supported')
    input = np.asarray(input)
    if np.iscomplexobj(input):
        raise TypeError('Complex type not supported')
    output = _get_output(output, input)
    if order in [0, 1]:
        output[...] = np.array(input)
    else:
        axis = _check_axis(axis, input.ndim)
        _nd_image.spline_filter1d(input, order, axis, output)
    return output    
    
def spline_filter(input, order=3, output=np.float64):
    if order < 2 or order > 5:
        raise RuntimeError('spline order not supported')
    input = np.asarray(input)
    if np.iscomplexobj(input):
        raise TypeError('Complex type not supported')
    output = _get_output(output, input)
    if order not in [0, 1] and input.ndim > 0:
        for axis in range(input.ndim):
            spline_filter1d(input, order, axis, output=output)
            input = output
    else:
        output[...] = input[...]
    return output

    
# --- extracted from:  scipy/ndimage/measurements.py ---    
    
def label(input, structure=None, output=None):
    input = np.asarray(input)
    if np.iscomplexobj(input):
        raise TypeError('Complex type not supported')
    if structure is None:
        structure = generate_binary_structure(input.ndim, 1)
    structure = np.asarray(structure, dtype=bool)
    if structure.ndim != input.ndim:
        raise RuntimeError('structure and input must have equal rank')
    for ii in structure.shape:
        if ii != 3:
            raise ValueError('structure dimensions must be equal to 3')

    # Use 32 bits if it's large enough for this image.
    # _ni_label.label()  needs two entries for background and
    # foreground tracking
    need_64bits = input.size >= (2**31 - 2)

    if isinstance(output, np.ndarray):
        if output.shape != input.shape:
            raise ValueError("output shape not correct")
        caller_provided_output = True
    else:
        caller_provided_output = False
        if output is None:
            output = np.empty(input.shape, np.intp if need_64bits else np.int32)
        else:
            output = np.empty(input.shape, output)

    # handle scalars, 0-dim arrays
    if input.ndim == 0 or input.size == 0:
        if input.ndim == 0:
            # scalar
            maxlabel = 1 if (input != 0) else 0
            output[...] = maxlabel
        else:
            # 0-dim
            maxlabel = 0
        if caller_provided_output:
            return maxlabel
        else:
            return output, maxlabel

    try:
        max_label = _ni_label._label(input, structure, output)
    except _ni_label.NeedMoreBits:
        # Make another attempt with enough bits, then try to cast to the
        # new type.
        tmp_output = np.empty(input.shape, np.intp if need_64bits else np.int32)
        max_label = _ni_label._label(input, structure, tmp_output)
        output[...] = tmp_output[...]
        if not np.all(output == tmp_output):
            # refuse to return bad results
            raise RuntimeError("insufficient bit-depth in requested output type")

    if caller_provided_output:
        # result was written in-place
        return max_label
    else:
        return output, max_label    

        
# --- extracted from:  scipy/ndimage/morphology.py ---           
        
def generate_binary_structure(rank, connectivity):
    if connectivity < 1:
        connectivity = 1
    if rank < 1:
        if connectivity < 1:
            return np.array(0, dtype=bool)
        else:
            return np.array(1, dtype=bool)
    output = np.fabs(np.indices([3] * rank) - 1)
    output = np.add.reduce(output, 0)
    return np.asarray(output <= connectivity, dtype=bool)


# --- extracted from:  scipy/_lib/_util.py ---

_MINPACK_LOCK = threading.RLock()

def _getargspec(func):
    # python 2.x 
    argspec = inspect.getargspec(func)
    if argspec.args[0] == 'self':
        argspec.args.pop(0)
    return argspec

    
# --- extracted from:  scipy/optimize/_lsq/least_squares.py ---    
    
def prepare_bounds(bounds, n):
    lb, ub = [np.asarray(b, dtype=float) for b in bounds]
    if lb.ndim == 0:
        lb = np.resize(lb, n)

    if ub.ndim == 0:
        ub = np.resize(ub, n)

    return lb, ub    
    
# --- extracted from:  scipy/optimize/minpack.py ---

def _check_func(checker, argname, thefunc, x0, args, numinputs,
                output_shape=None):
    res = np.atleast_1d(thefunc(*((x0[:numinputs],) + args)))
    if (output_shape is not None) and (shape(res) != output_shape):
        if (output_shape[0] != 1):
            if len(output_shape) > 1:
                if output_shape[1] == 1:
                    return shape(res)
            msg = "%s: there is a mismatch between the input and output " \
                  "shape of the '%s' argument" % (checker, argname)
            func_name = getattr(thefunc, '__name__', None)
            if func_name:
                msg += " '%s'." % func_name
            else:
                msg += "."
            msg += 'Shape should be %s but it is %s.' % (output_shape, shape(res))
            raise TypeError(msg)
    if np.issubdtype(res.dtype, np.inexact):
        dt = res.dtype
    else:
        dt = dtype(float)
    return np.shape(res), dt
    
def leastsq(func, x0, args=(), Dfun=None, full_output=0,
            col_deriv=0, ftol=1.49012e-8, xtol=1.49012e-8,
            gtol=0.0, maxfev=0, epsfcn=None, factor=100, diag=None):

    x0 = np.asarray(x0).flatten()
    n = len(x0)
    if not isinstance(args, tuple):
        args = (args,)
    shape, dtype = _check_func('leastsq', 'func', func, x0, args, n)
    m = shape[0]
    if n > m:
        raise TypeError('Improper input: N=%s must not exceed M=%s' % (n, m))
    if epsfcn is None:
        epsfcn = np.finfo(dtype).eps
    if Dfun is None:
        if maxfev == 0:
            maxfev = 200*(n + 1)
        with _MINPACK_LOCK:
            retval = _minpack._lmdif(func, x0, args, full_output, ftol, xtol,
                                     gtol, maxfev, epsfcn, factor, diag)
    else:
        if col_deriv:
            _check_func('leastsq', 'Dfun', Dfun, x0, args, n, (n, m))
        else:
            _check_func('leastsq', 'Dfun', Dfun, x0, args, n, (m, n))
        if maxfev == 0:
            maxfev = 100 * (n + 1)
        with _MINPACK_LOCK:
            retval = _minpack._lmder(func, Dfun, x0, args, full_output,
                                     col_deriv, ftol, xtol, gtol, maxfev,
                                     factor, diag)

    errors = {0: ["Improper input parameters.", TypeError],
              1: ["Both actual and predicted relative reductions "
                  "in the sum of squares\n  are at most %f" % ftol, None],
              2: ["The relative error between two consecutive "
                  "iterates is at most %f" % xtol, None],
              3: ["Both actual and predicted relative reductions in "
                  "the sum of squares\n  are at most %f and the "
                  "relative error between two consecutive "
                  "iterates is at \n  most %f" % (ftol, xtol), None],
              4: ["The cosine of the angle between func(x) and any "
                  "column of the\n  Jacobian is at most %f in "
                  "absolute value" % gtol, None],
              5: ["Number of calls to function has reached "
                  "maxfev = %d." % maxfev, ValueError],
              6: ["ftol=%f is too small, no further reduction "
                  "in the sum of squares\n  is possible.""" % ftol,
                  ValueError],
              7: ["xtol=%f is too small, no further improvement in "
                  "the approximate\n  solution is possible." % xtol,
                  ValueError],
              8: ["gtol=%f is too small, func(x) is orthogonal to the "
                  "columns of\n  the Jacobian to machine "
                  "precision." % gtol, ValueError],
              'unknown': ["Unknown error.", TypeError]}

    info = retval[-1]    # The FORTRAN return value

    if info not in [1, 2, 3, 4] and not full_output:
            try:
                raise errors[info][1](errors[info][0])
            except KeyError:
                raise errors['unknown'][1](errors['unknown'][0])

    mesg = errors[info][0]
    if full_output:
        cov_x = None
        if info in [1, 2, 3, 4]:
            from numpy.dual import inv
            perm = np.take(np.eye(n), retval[1]['ipvt'] - 1, 0)
            r = np.triu(np.transpose(retval[1]['fjac'])[:n, :])
            R = np.dot(r, perm)
            try:
                cov_x = inv(dot(transpose(R), R))
            except:
                pass
        return (retval[0], cov_x) + retval[1:-1] + (mesg, info)
    else:
        return (retval[0], info)

def _wrap_func(func, xdata, ydata, transform):
    if transform is None:
        def func_wrapped(params):
            return func(xdata, *params) - ydata
    elif transform.ndim == 1:
        def func_wrapped(params):
            return transform * (func(xdata, *params) - ydata)
    else:
        def func_wrapped(params):
            return solve_triangular(transform, func(xdata, *params) - ydata, lower=True)
    return func_wrapped

def curve_fit(f, xdata, ydata, p0=None, sigma=None, absolute_sigma=False,
              check_finite=True, bounds=(-np.inf, np.inf), method=None,
              jac=None, **kwargs):
              
    if p0 is None:
        # determine number of parameters by inspecting the function
        #from scipy._lib._util import getargspec_no_self as _getargspec
        args, varargs, varkw, defaults = _getargspec(f)
        if len(args) < 2:
            raise ValueError("Unable to determine number of fit parameters.")
        n = len(args) - 1
    else:
        p0 = np.atleast_1d(p0)
        n = p0.size

    lb, ub = prepare_bounds(bounds, n)
    if p0 is None:
        p0 = _initialize_feasible(lb, ub)

    bounded_problem = np.any((lb > -np.inf) | (ub < np.inf))
    if method is None:
        if bounded_problem:
            method = 'trf'
        else:
            method = 'lm'

    if method == 'lm' and bounded_problem:
        raise ValueError("Method 'lm' only works for unconstrained problems. "
                         "Use 'trf' or 'dogbox' instead.")

    # NaNs can not be handled
    if check_finite:
        ydata = np.asarray_chkfinite(ydata)
    else:
        ydata = np.asarray(ydata)

    if isinstance(xdata, (list, tuple, np.ndarray)):
        # `xdata` is passed straight to the user-defined `f`, so allow
        # non-array_like `xdata`.
        if check_finite:
            xdata = np.asarray_chkfinite(xdata)
        else:
            xdata = np.asarray(xdata)

    # Determine type of sigma
    if sigma is not None:
        sigma = np.asarray(sigma)

        # if 1-d, sigma are errors, define transform = 1/sigma
        if sigma.shape == (ydata.size, ):
            transform = 1.0 / sigma
        # if 2-d, sigma is the covariance matrix,
        # define transform = L such that L L^T = C
        elif sigma.shape == (ydata.size, ydata.size):
            try:
                # scipy.linalg.cholesky requires lower=True to return L L^T = A
                transform = cholesky(sigma, lower=True)
            except:
                raise ValueError("`sigma` must be positive definite.")
        else:
            raise ValueError("`sigma` has incorrect shape.")
    else:
        transform = None

    func = _wrap_func(f, xdata, ydata, transform)
    if callable(jac):
        jac = _wrap_jac(jac, xdata, transform)
    elif jac is None and method != 'lm':
        jac = '2-point'

    if method == 'lm':
        # Remove full_output from kwargs, otherwise we're passing it in twice.
        return_full = kwargs.pop('full_output', False)
        res = leastsq(func, p0, Dfun=jac, full_output=1, **kwargs)
        popt, pcov, infodict, errmsg, ier = res
        cost = np.sum(infodict['fvec'] ** 2)
        if ier not in [1, 2, 3, 4]:
            raise RuntimeError("Optimal parameters not found: " + errmsg)
    else:
        # Rename maxfev (leastsq) to max_nfev (least_squares), if specified.
        if 'max_nfev' not in kwargs:
            kwargs['max_nfev'] = kwargs.pop('maxfev', None)

        res = least_squares(func, p0, jac=jac, bounds=bounds, method=method,
                            **kwargs)

        if not res.success:
            raise RuntimeError("Optimal parameters not found: " + res.message)

        cost = 2 * res.cost  # res.cost is half sum of squares!
        popt = res.x

        # Do Moore-Penrose inverse discarding zero singular values.
        _, s, VT = svd(res.jac, full_matrices=False)
        threshold = np.finfo(float).eps * max(res.jac.shape) * s[0]
        s = s[s > threshold]
        VT = VT[:s.size]
        pcov = np.dot(VT.T / s**2, VT)
        return_full = False

    warn_cov = False
    if pcov is None:
        # indeterminate covariance
        pcov = np.zeros((len(popt), len(popt)), dtype=float)
        pcov.fill(np.inf)
        warn_cov = True
    elif not absolute_sigma:
        if ydata.size > p0.size:
            s_sq = cost / (ydata.size - p0.size)
            pcov = pcov * s_sq
        else:
            pcov.fill(np.inf)
            warn_cov = True

    if return_full:
        return popt, pcov, infodict, errmsg, ier
    else:
        return popt, pcov



    
                              
# ----------------- test --------------------------- 
'''
test = np.random.random_sample ((28,28))
test = zoom(test[:,:],[0.25,0.25],order=1) # downsample
print (test)
print ('\n')      
test = median_filter(test, size = (3,3)) # median filter
s = 2; w = 4; t = (((w - 1)/2)-0.5)/s
test[:,:] = gaussian_filter(test[:,:], sigma=s, truncate=t)
print (test)
print ('\n')
print ('\n') 
mask = test > 0.5; mask = mask.astype(np.int16)
print (mask)
print ('\n') 
s = [[1,1,1],[1,1,1],[1,1,1]]
labeled_mask, num_clusters = label(mask, structure=s)
print (labeled_mask)
'''       