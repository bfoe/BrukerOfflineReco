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
# to be copied to the directory of this script
#

import numpy
import sys
import math
import _nd_image
import _ni_label

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
        output = numpy.zeros(shape, dtype=input.dtype.name)
    elif type(output) in [type(type), type(numpy.zeros((4,)).dtype)]:
        output = numpy.zeros(shape, dtype=output)
    elif type(output) in string_types:
        output = numpy.typeDict[output]
        output = numpy.zeros(shape, dtype=output)
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
    input = numpy.asarray(input)
    if numpy.iscomplexobj(input):
        raise TypeError('Complex type not supported')
    origins = _normalize_sequence(origin, input.ndim)
    if footprint is None:
        if size is None:
            raise RuntimeError("no footprint or filter size provided")
        sizes = _normalize_sequence(size, input.ndim)
        footprint = numpy.ones(sizes, dtype=bool)
    else:
        footprint = numpy.asarray(footprint, dtype=bool)
    fshape = [ii for ii in footprint.shape if ii > 0]
    if len(fshape) != input.ndim:
        raise RuntimeError('filter footprint array has incorrect shape.')
    for origin, lenf in zip(origins, fshape):
        if (lenf // 2 + origin < 0) or (lenf // 2 + origin >= lenf):
            raise ValueError('invalid origin')
    if not footprint.flags.contiguous:
        footprint = footprint.copy()
    filter_size = numpy.where(footprint, 1, 0).sum()
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
        output, return_value = _get_output(output, input)
        mode = _extend_mode_to_code(mode)
        _nd_image.rank_filter(input, rank, footprint, output, mode, cval,
                              origins)
        return return_value

        
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
    input = numpy.asarray(input)
    if numpy.iscomplexobj(input):
        raise TypeError('Complex type not supported')
    output, return_value = _get_output(output, input)
    weights = numpy.asarray(weights, dtype=numpy.float64)
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
    return return_value            
        
def median_filter(input, size=None, footprint=None, output=None,
                  mode="reflect", cval=0.0, origin=0):  
    return _rank_filter(input, 0, size, footprint, output, mode, cval,
                        origin, 'median')  

def gaussian_filter(input, sigma, order=0, output=None,
      mode="reflect", cval=0.0, truncate=4.0):
    input = numpy.asarray(input)
    output, return_value = _get_output(output, input)
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
    return return_value


# --- extracted from:  scipy/ndimage/interpolation.py ---
    
def zoom(input, zoom, output=None, order=3, mode='constant', cval=0.0,
         prefilter=True):
    if order < 0 or order > 5:
        raise RuntimeError('spline order not supported')
    input = numpy.asarray(input)
    if numpy.iscomplexobj(input):
        raise TypeError('Complex type not supported')
    if input.ndim < 1:
        raise RuntimeError('input and output rank must be > 0')
    mode = _extend_mode_to_code(mode)
    if prefilter and order > 1:
        filtered = spline_filter(input, order, output=numpy.float64)
    else:
        filtered = input
    zoom = _normalize_sequence(zoom, input.ndim)
    output_shape = tuple(
            [int(round(ii * jj)) for ii, jj in zip(input.shape, zoom)])

    output_shape_old = tuple(
            [int(ii * jj) for ii, jj in zip(input.shape, zoom)])
    if output_shape != output_shape_old:
        warnings.warn(
                "From scipy 0.13.0, the output shape of zoom() is calculated "
                "with round() instead of int() - for these inputs the size of "
                "the returned array has changed.", UserWarning)

    zoom_div = numpy.array(output_shape, float) - 1
    # Zooming to infinite values is unpredictable, so just choose
    # zoom factor 1 instead
    zoom = numpy.divide(numpy.array(input.shape) - 1, zoom_div,
                        out=numpy.ones_like(input.shape, dtype=numpy.float64),
                        where=zoom_div != 0)

    output, return_value = _get_output(output, input,
                                                   shape=output_shape)
    zoom = numpy.ascontiguousarray(zoom)
    _nd_image.zoom_shift(filtered, zoom, None, output, order, mode, cval)
    return return_value

def rotate(input, angle, axes=(1, 0), reshape=True, output=None, order=3,
           mode='constant', cval=0.0, prefilter=True):
    input = numpy.asarray(input)
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
    angle = numpy.pi / 180 * angle
    m11 = math.cos(angle)
    m12 = math.sin(angle)
    m21 = -math.sin(angle)
    m22 = math.cos(angle)
    matrix = numpy.array([[m11, m12],
                          [m21, m22]], dtype=numpy.float64)
    iy = input.shape[axes[0]]
    ix = input.shape[axes[1]]
    if reshape:
        mtrx = numpy.array([[m11, -m21],
                            [-m12, m22]], dtype=numpy.float64)
        minc = [0, 0]
        maxc = [0, 0]
        coor = numpy.dot(mtrx, [0, ix])
        minc, maxc = _minmax(coor, minc, maxc)
        coor = numpy.dot(mtrx, [iy, 0])
        minc, maxc = _minmax(coor, minc, maxc)
        coor = numpy.dot(mtrx, [iy, ix])
        minc, maxc = _minmax(coor, minc, maxc)
        oy = int(maxc[0] - minc[0] + 0.5)
        ox = int(maxc[1] - minc[1] + 0.5)
    else:
        oy = input.shape[axes[0]]
        ox = input.shape[axes[1]]
    offset = numpy.zeros((2,), dtype=numpy.float64)
    offset[0] = float(oy) / 2.0 - 0.5
    offset[1] = float(ox) / 2.0 - 0.5
    offset = numpy.dot(matrix, offset)
    tmp = numpy.zeros((2,), dtype=numpy.float64)
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
        size = numpy.product(input.shape, axis=0)
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
    input = numpy.asarray(input)
    if numpy.iscomplexobj(input):
        raise TypeError('Complex type not supported')
    if output_shape is None:
        output_shape = input.shape
    if input.ndim < 1 or len(output_shape) < 1:
        raise RuntimeError('input and output rank must be > 0')
    mode = _extend_mode_to_code(mode)
    if prefilter and order > 1:
        filtered = spline_filter(input, order, output=numpy.float64)
    else:
        filtered = input
    output = _get_output(output, input,
                                                   shape=output_shape)
    matrix = numpy.asarray(matrix, dtype=numpy.float64)
    if matrix.ndim not in [1, 2] or matrix.shape[0] < 1:
        raise RuntimeError('no proper affine matrix provided')
    if (matrix.ndim == 2 and matrix.shape[1] == input.ndim + 1 and
            (matrix.shape[0] in [input.ndim, input.ndim + 1])):
        if matrix.shape[0] == input.ndim + 1:
            exptd = [0] * input.ndim + [1]
            if not numpy.all(matrix[input.ndim] == exptd):
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
    offset = numpy.asarray(offset, dtype=numpy.float64)
    if offset.ndim != 1 or offset.shape[0] < 1:
        raise RuntimeError('no proper offset provided')
    if not offset.flags.contiguous:
        offset = offset.copy()
    if matrix.ndim == 1:
        warnings.warn(
            "The behaviour of affine_transform with a one-dimensional "
            "array supplied for the matrix parameter has changed in "
            "scipy 0.18.0."
        )
        _nd_image.zoom_shift(filtered, matrix, offset/matrix, output, order,
                             mode, cval)
    else:
        _nd_image.geometric_transform(filtered, None, None, matrix, offset,
                                      output, order, mode, cval, None, None)
    return output

def spline_filter1d(input, order=3, axis=-1, output=numpy.float64):
    if order < 0 or order > 5:
        raise RuntimeError('spline order not supported')
    input = numpy.asarray(input)
    if numpy.iscomplexobj(input):
        raise TypeError('Complex type not supported')
    output = _get_output(output, input)
    if order in [0, 1]:
        output[...] = numpy.array(input)
    else:
        axis = _check_axis(axis, input.ndim)
        _nd_image.spline_filter1d(input, order, axis, output)
    return output    
    
def spline_filter(input, order=3, output=numpy.float64):
    if order < 2 or order > 5:
        raise RuntimeError('spline order not supported')
    input = numpy.asarray(input)
    if numpy.iscomplexobj(input):
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
    input = numpy.asarray(input)
    if numpy.iscomplexobj(input):
        raise TypeError('Complex type not supported')
    if structure is None:
        structure = generate_binary_structure(input.ndim, 1)
    structure = numpy.asarray(structure, dtype=bool)
    if structure.ndim != input.ndim:
        raise RuntimeError('structure and input must have equal rank')
    for ii in structure.shape:
        if ii != 3:
            raise ValueError('structure dimensions must be equal to 3')

    # Use 32 bits if it's large enough for this image.
    # _ni_label.label()  needs two entries for background and
    # foreground tracking
    need_64bits = input.size >= (2**31 - 2)

    if isinstance(output, numpy.ndarray):
        if output.shape != input.shape:
            raise ValueError("output shape not correct")
        caller_provided_output = True
    else:
        caller_provided_output = False
        if output is None:
            output = numpy.empty(input.shape, numpy.intp if need_64bits else numpy.int32)
        else:
            output = numpy.empty(input.shape, output)

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
        tmp_output = numpy.empty(input.shape, numpy.intp if need_64bits else numpy.int32)
        max_label = _ni_label._label(input, structure, tmp_output)
        output[...] = tmp_output[...]
        if not numpy.all(output == tmp_output):
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
            return numpy.array(0, dtype=bool)
        else:
            return numpy.array(1, dtype=bool)
    output = numpy.fabs(numpy.indices([3] * rank) - 1)
    output = numpy.add.reduce(output, 0)
    return numpy.asarray(output <= connectivity, dtype=bool)
        
                              
# ----------------- test --------------------------- 
'''
test = numpy.random.random_sample ((28,28))
test = zoom(test[:,:],[0.25,0.25],order=1) # downsample
print (test)
print ('\n')      
test = median_filter(test, size = (3,3)) # median filter
s = 2; w = 4; t = (((w - 1)/2)-0.5)/s
test[:,:] = gaussian_filter(test[:,:], sigma=s, truncate=t)
print (test)

print ('\n')
print ('\n') 
mask = test > 0.5; mask = mask.astype(numpy.int16)
print (mask)
print ('\n') 
s = [[1,1,1],[1,1,1],[1,1,1]]
labeled_mask, num_clusters = label(mask, structure=s)
print (labeled_mask)
'''       