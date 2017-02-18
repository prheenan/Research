#! /usr/bin/python
# -*- coding: utf-8 -*-
'''
This module contains the second stage posprocessing tools.
'''
__author__  = "Charles Roduit <charles.roduit@gmail.com>"
__date__ = "28.12.2008"
__license__ = "GNU Public License (GPL) version 3"
__version__ = "0.1"

from warnings import warn

import numpy
from scipy import stats

from misc import ismasked

def get_array_distance(source_array, distance):
    r'''
    Create a 3D array that contains the array series of the pixels distances 
    from each marked pixels.
    
    for a distance of 1, the source and result arrays look like this :
    
    ::
    
       source_matrix   result[:,:,0]   result[:, :, 1]
         0 0 0 0 0       0 0 0 0 0       0 0 0 0 0   
         0 0 1 0 0       0 1 0 1 0       0 0 0 0 0   
         0 0 0 0 0       0 0 0 0 0       0 0 0 0 0   
         0 1 0 0 0       0 0 0 0 0       1 0 1 0 0   
         0 0 0 0 0       0 0 0 0 0       0 0 0 0 0   
    '''
    
    [pos_x, pos_y] = source_array.nonzero()
    all_drops = drop_seed(source_array, distance)
    all_single_drops = numpy.empty((source_array.shape[0], 
                                   source_array.shape[1], 
                                   len(pos_x)))
    all_single_seeds = numpy.empty_like(all_single_drops)

    for item in range(len(pos_x)):
        this_seed = numpy.zeros(source_array.shape)
        this_seed[pos_x[item], pos_y[item]] = 1
        this_drop = drop_seed(this_seed, distance)
        all_single_drops[:, :, item] = this_drop * all_drops
        all_single_seeds[:, :, item] = this_seed
    return [all_single_seeds, all_single_drops]
    
def drop_seed(seed_matrix, size):
    '''
    This function returns a similar matrix as dist_seed, except that it takes
    care of the overlapping
    '''
    
    drop_matrix = numpy.zeros_like(seed_matrix)
    for dist in range(1, size):
        drop_matrix = drop_matrix + dist_seed(seed_matrix, dist, where = 'xy')
        
    negative_matrix = numpy.logical_not(drop_matrix + seed_matrix)
    drop_matrix = drop_matrix + dist_seed(seed_matrix, size)
    border_matrix = drop_matrix * negative_matrix
    
    return border_matrix
    
def dist_seed(seed_matrix, distance, where = 'x'):
    '''
    This function returns a daughter matrix that contains the position of pixels
    at the given distance from the non-zero pixel of the source matrix.
    
    ::
    
         Source            Result
       0 0 0 0 0         0 0 0 0 0 
       0 0 0 0 0         0 0 0 0 0 
       0 0 1 0 0   ==>   0 1 0 1 0 
       0 0 0 0 0         0 0 0 0 0 
       0 0 0 0 0         0 0 0 0 0 
    
    '''
    
    [pos_x, pos_y] = seed_matrix.nonzero()
    dist_matrix = numpy.zeros_like(seed_matrix)
    if where == 'x':
        to_left = pos_x - distance
        to_right = pos_x + distance
        # test if the pixels are in or outside the matrix
        limit_left = ( to_left >= 0 ) * ( to_left < len(seed_matrix) )
        limit_right = ( to_right >= 0 ) * ( to_right < len(seed_matrix) )
        # Take the pixels that are in the matrix
        pos_x_left = to_left[limit_left]
        pos_x_right = to_right[limit_right]
        pos_y_left = pos_y[limit_left]
        pos_y_right = pos_y[limit_right]        
        
        # Generate the matrix
        dist_matrix[pos_x_left, pos_y_left] = 1
        dist_matrix[pos_x_right, pos_y_right] = 1
    elif where == 'xy':
        for elm_nbr in range(len(pos_x)):
            this_x = pos_x[elm_nbr] + distance
            this_y = pos_y[elm_nbr]
            if this_x < dist_matrix.shape[0] and this_y < dist_matrix.shape[1]:
                dist_matrix[this_x, this_y] = 1
            # from 0 to pi/2
            while this_y < pos_y[elm_nbr] + distance:
                if (this_x < dist_matrix.shape[0] and this_x >= 0 and
                   this_y < dist_matrix.shape[1] and this_y >= 0):
                    dist_matrix[this_x, this_y] = 1
                this_y = this_y + 1
                if numpy.sqrt((this_x - pos_x[elm_nbr]) ** 2 + 
                              (this_y - pos_y[elm_nbr]) ** 2) > distance :
                    this_x = this_x - 1
                if numpy.sqrt((this_x - pos_x[elm_nbr]) ** 2 + 
                              (this_y - pos_y[elm_nbr]) ** 2) > distance :
                    this_y = this_y - 1
            # from pi/2 to pi
            while this_x > pos_x[elm_nbr] - distance:
                if (this_x < dist_matrix.shape[0] and this_x >= 0 and
                   this_y < dist_matrix.shape[1] and this_y >= 0):
                    dist_matrix[this_x, this_y] = 1
                this_x = this_x - 1
                if numpy.sqrt((this_x - pos_x[elm_nbr]) ** 2 + 
                              (this_y - pos_y[elm_nbr]) ** 2) > distance :
                    this_y = this_y - 1
                if numpy.sqrt((this_x - pos_x[elm_nbr]) ** 2 + 
                              (this_y - pos_y[elm_nbr]) ** 2) > distance :
                    this_x = this_x + 1
            # from pi to 3pi/2
            while this_y > pos_y[elm_nbr] - distance:
                if (this_x < dist_matrix.shape[0] and this_x >= 0 and
                   this_y < dist_matrix.shape[1] and this_y >= 0):
                    dist_matrix[this_x, this_y] = 1
                this_y = this_y - 1
                if numpy.sqrt((this_x - pos_x[elm_nbr]) ** 2 + 
                              (this_y - pos_y[elm_nbr]) ** 2) > distance :
                    this_x = this_x + 1
                if numpy.sqrt((this_x - pos_x[elm_nbr]) ** 2 + 
                              (this_y - pos_y[elm_nbr]) ** 2) > distance :
                    this_y = this_y + 1
            # from 3pi/2 to 2pi
            while this_x < pos_x[elm_nbr] + distance:
                if (this_x < dist_matrix.shape[0] and this_x >= 0 and
                   this_y < dist_matrix.shape[1] and this_y >= 0):
                    dist_matrix[this_x, this_y] = 1
                this_x = this_x + 1
                if numpy.sqrt((this_x - pos_x[elm_nbr]) ** 2 + 
                              (this_y - pos_y[elm_nbr]) ** 2) > distance :
                    this_y = this_y + 1
                if numpy.sqrt((this_x - pos_x[elm_nbr]) ** 2 + 
                              (this_y - pos_y[elm_nbr]) ** 2) > distance :
                    this_x = this_x - 1
            
            #while this_y 
            pass
    return dist_matrix
 
def rel_values(value_array, dot_array, max_dist = 3, coord=False):
    '''
    This function computes the relative values at some points defined by the
    dot_array.
    
    dot_array = n * m array that contains the position where to compute relative
                values.
                
    value_array = n * m array that contains the values.
    
    max_dist = integer, to which distance to compute (= 3 by default)
    
    Output :
    
    result[dot_nb, rel_val] = result[0, 1] is the absolute value of the first
                                           dot
                                           
                              result[1:max_dist, 1] is the relative values of
                                                    the first dot
    '''
    if type(value_array) == numpy.ma.core.MaskedArray:
        #print type(dot_array)
        seed_array = dot_array * (1 - value_array.mask)
    else:
        seed_array = dot_array
    if coord:
        # we record the coordinate of the events.
        coord_val = []
    # we compute the distance matrices... 
    # First val is the absolute value.
    result = numpy.empty((max_dist + 1, seed_array.astype(bool).sum()))
    ret_array = []
    for distance in range(max_dist + 1):
        [single_dot, rel_array] = get_array_distance(seed_array, distance)
        for dot_nb in range(single_dot.shape[2]):
            # we compute the value at the dot position...
            mask_dot = single_dot[:, :, dot_nb].astype(bool)
            dot_val = value_array[mask_dot].sum()
            #print type(dot_val)
            #print "distance %i, dot number %i, dot val %f"%(distance, dot_nb, dot_val)
            if not distance:
                # we want the first value to be the value at the position
                result[0, dot_nb] = dot_val
                ret_array.append(single_dot)
                if coord:
                    coord_val.append([mask_dot.nonzero()[0][0], 
                                      mask_dot.nonzero()[1][0]])
            else:
                # ... and at the distant positions ...
                mask_rel = rel_array[:, :, dot_nb].astype(bool)
                dist_val = value_array[mask_rel].sum()
                # ... with the number of distant positions ...
                try:
                    count_dist = len(value_array[mask_rel])
                    # ... and finally the relative value :
                    # First exclude the masked and negative datas
                    if ismasked(dist_val) and dist_val.mask:
                        rel_value = numpy.nan
                    elif dot_val <= 0 or dist_val <= 0:
                        rel_value = numpy.nan
                    else:
                        # to finally compute the relative values.
                        rel_value = numpy.log(dot_val / dist_val * count_dist)
                except ZeroDivisionError:
                    rel_value = numpy.inf
                if rel_value == numpy.inf:
                    result[distance, dot_nb] = numpy.nan
                else:
                    result[distance, dot_nb] = rel_value
                ret_array.append(rel_array)
    # we mask the nan
    result = numpy.ma.masked_array(result, mask = numpy.isnan(result))
    if coord:
        return result, coord_val
    else:
        return result#, result_summarize
def add_stiffness_to_event(event_list, rel_stiff, coordinate):
    for event in event_list:
        try:
            __index = coordinate.index([event['pos_x'], event['pos_y']])
        except ValueError:
            event['Stiffness'] = numpy.ma.array([], mask=[])
        else:
            event['Stiffness'] = rel_stiff[:, __index]
    return event_list
def rand_dot(dot_array, nb_dot = None):
    '''
    Creates an array with dots. The dots are placed randomly, with the only
    constrain to place in a free position, given by the input dot array.
    '''
    nb_init_dot = len(dot_array.nonzero()[0])
    
    if nb_dot == None:
        nb_dot = nb_init_dot
    if nb_dot > dot_array.size - nb_init_dot:
        warn('Not enough place to generate random events.', RuntimeWarning)
        # have to reach an error...
        return numpy.zeros_like(dot_array)
    result_array = numpy.zeros(dot_array.shape)
    while result_array.sum() < nb_dot:
        pos_x = numpy.floor(dot_array.shape[0] * numpy.random.rand())
        pos_y = numpy.floor(dot_array.shape[1] * numpy.random.rand())
        if not(dot_array[pos_x, pos_y] or result_array[pos_x, pos_y]):
            result_array[pos_x, pos_y] = 1
    return result_array

def resize_array(old_array, new_shape):
    '''
    Create a new array with the given size and with the data in old array.
    '''
    # The new arrays are made of NaN to differentiate data from empty blocks
    new_array = numpy.zeros(new_shape, numpy.float) * numpy.nan
    # We reposition the points in the new array.
    for pos_y in range(old_array.shape[1]):
        for pos_x in range(old_array.shape[0]):
            curve_nb = pos_x + pos_y * old_array.shape[0]
            n_pos_x = (curve_nb) % new_shape[0]
            n_pos_y = (curve_nb) / new_shape[0]
            new_array[n_pos_x, n_pos_y] = old_array[pos_x, pos_y]
    return new_array

def kde_fit(data, nbr_div = 128):
    '''
    Automatically detects best gaussian fits for the data.
    
    For example, let's create a data set of to normal random values :
    
    >>> import pylab
    
    >>> data1 = numpy.random.randn(100)
    >>> data = numpy.zeros(200)
    >>> data[:100] = data1+1
    >>> data[100:] = data1+5
    
    Now, let's compute the kde fit :
    
    >>> [fit, maxima, minima] = kde_fit(data)
    
    And we can plot the result :
    
    >>> [x, y, patch] = pylab.hist(data)
    >>> conversion_factor = max(x)/max(fit[1])
    >>> figure = pylab.figure()
    >>> axis = figure.add_subplot(111)
    >>> pl = pylab.plot(fit[0], fit[1] * conversion_factor)
    >>> axis.hold(True)
    >>> pl = pylab.plot(maxima[0], maxima[1] * conversion_factor, 'go')
    >>> pl = pylab.plot(minima[0], minima[1] * conversion_factor, 'ro')
    >>> figure.savefig('figure/kde_fit.png')
    
    The resulting graph will be an hisogram of the data set with the kde fit
    plotted in foreground. Two gree dots represents the maximas and the red dot
    represent the minima.
    '''
    kde = stats.kde.gaussian_kde(data)
    min_data = min(data)
    max_data = max(data)
    nbr_step = (max_data - min_data) / nbr_div
    coord = numpy.arange(min_data, max_data, nbr_step)
    gaussian_fit = kde(coord)
    # After the gaussian fit, we want to finds the maximums
    deriv = numpy.diff(gaussian_fit)
    #print deriv
    extrema = (abs(deriv) == deriv) * 1
    extrema = numpy.diff(extrema)
    #print extremums
    maxima = numpy.nonzero(extrema == -1)[0] + 1 # add one to correspond to
    minima = numpy.nonzero(extrema == 1)[0] + 1  # the maximum
    return [[coord, gaussian_fit], 
            [coord[maxima], numpy.ones(len(maxima)) * gaussian_fit[maxima]],
            [coord[minima], numpy.ones(len(minima)) * gaussian_fit[minima]]]
def compute_gauss(data, minima, maxima):
    '''
    compute the characteristic of the distribution : mean, std, and number of
    elements.
    '''
    mean, std, count = [], [], []
    position = []
    if not len(minima) or minima[0] > maxima[0]:
        position.append(min(data))
    for minimum in minima:
        position.append(minimum)
    if not len(minima) or minima[-1] < maxima[-1]:
        position.append(max(data))
    for begin, end in zip(position[:-1], position[1:]):
        index = (data >= begin) * (data < end)
        if len(data[index]):
            try:
                mean.append(numpy.mean(data[index]))
                std.append(numpy.std(data[index]))
                count.append(len(data[index]))
            except TypeError:
                # In windows, numpy.nanmax seems to not work in this case...
                index = index * numpy.array((numpy.isnan(data)-1), dtype=bool)
                mean.append(numpy.mean(data[index]))
                std.append(numpy.std(data[index]))
                count.append(len(data[index]))
        else:
            # Array is empty...
            mean.append(numpy.nan)
            std.append(numpy.nan)
            count.append(len(data[index]))
    return [mean, std, count]
def find_gauss_fit(data, nbr_div = 128):
    '''
    Automatically detects best gaussian fits for the data, and retun the
    gaussian properties of each detected peak.
    
    For example, let's create a data set of to normal random values :
    
    >>> import pylab
    
    >>> data1 = numpy.random.randn(100)
    >>> data = numpy.zeros(200)
    >>> data[:100] = data1+1
    >>> data[100:] = data1+5
    
    Now, we can compute the gauss fit of these data :
    
    >>> fit = find_gauss_fit(data)
    
    And plot the result :
    
    >>> [x, y, patch] = pylab.hist(data)
    >>> conversion_factor = max(x)/max(fit['kde'][1])
    >>> figure = pylab.figure()
    >>> axis = figure.add_subplot(111)
    >>> pl = axis.hist(data)
    >>> axis.hold(True)
    >>> pl = axis.plot(fit['kde'][0], fit['kde'][1] * conversion_factor)
    >>> pl = axis.plot(fit['maxima'][0], fit['maxima'][1] * conversion_factor, 'go')
    >>> pl = axis.plot(fit['mean'], fit['maxima'][1] * conversion_factor, 'ro')
    >>> pl = axis.errorbar(fit['mean'], fit['maxima'][1] * conversion_factor, xerr=fit['std'], fmt=None, ecolor='red')
    >>> figure.savefig('figure/gauss_fit.png')
    
    The resulting plot shows the data as an histogram with the kde fit in the
    foreground.
    The detected maxima are shown as two green dots and the mean as red dot with
    the horizontal error bar describing the std deviation.
    '''
    data = data[numpy.isfinite(data)]
    if len(data):
        [fit, maxima, minima] = kde_fit(data, nbr_div)
        [mean_lst, std_lst, count_lst] = compute_gauss(data, minima[0], maxima[0])
    else:
        [fit, maxima, minima, mean_lst, std_lst, count_lst] = [[] for i in range(6)]
           
    return {'kde' : fit, 'mean' : mean_lst, 'std' : std_lst, 'count' :
            count_lst, 'maxima' : maxima}

def tomography_array(stiffness, floor=None):
    """
        Generate the stiffness tomography array.
        
        >>> raw_array = numpy.mgrid[0:5, 0:5, 0:5][2]
        >>> floor = numpy.mgrid[0:5, 0:5, 0:5][1][0]
        >>> tomo = tomography_array(raw_array, floor)
        >>> print tomo
        [[[4.0 3.0 2.0 1.0 0.0 -- -- -- --]
          [4.0 3.0 2.0 1.0 0.0 -- -- -- --]
          [4.0 3.0 2.0 1.0 0.0 -- -- -- --]
          [4.0 3.0 2.0 1.0 0.0 -- -- -- --]
          [4.0 3.0 2.0 1.0 0.0 -- -- -- --]]
        <BLANKLINE>
         [[-- 4.0 3.0 2.0 1.0 0.0 -- -- --]
          [-- 4.0 3.0 2.0 1.0 0.0 -- -- --]
          [-- 4.0 3.0 2.0 1.0 0.0 -- -- --]
          [-- 4.0 3.0 2.0 1.0 0.0 -- -- --]
          [-- 4.0 3.0 2.0 1.0 0.0 -- -- --]]
        <BLANKLINE>
         [[-- -- 4.0 3.0 2.0 1.0 0.0 -- --]
          [-- -- 4.0 3.0 2.0 1.0 0.0 -- --]
          [-- -- 4.0 3.0 2.0 1.0 0.0 -- --]
          [-- -- 4.0 3.0 2.0 1.0 0.0 -- --]
          [-- -- 4.0 3.0 2.0 1.0 0.0 -- --]]
        <BLANKLINE>
         [[-- -- -- 4.0 3.0 2.0 1.0 0.0 --]
          [-- -- -- 4.0 3.0 2.0 1.0 0.0 --]
          [-- -- -- 4.0 3.0 2.0 1.0 0.0 --]
          [-- -- -- 4.0 3.0 2.0 1.0 0.0 --]
          [-- -- -- 4.0 3.0 2.0 1.0 0.0 --]]
        <BLANKLINE>
         [[-- -- -- -- 4.0 3.0 2.0 1.0 0.0]
          [-- -- -- -- 4.0 3.0 2.0 1.0 0.0]
          [-- -- -- -- 4.0 3.0 2.0 1.0 0.0]
          [-- -- -- -- 4.0 3.0 2.0 1.0 0.0]
          [-- -- -- -- 4.0 3.0 2.0 1.0 0.0]]]
    """
    # The piezo array contain the height of the tip at the end of the
    # indentation. It correspond then to the last slice of the stiffness
    # array.
    if floor is None:
        floor = numpy.zeros(stiffness.shape[0:2])
    max_floor = floor.max()
    
    new_array = numpy.zeros((stiffness.shape[0], stiffness.shape[1], 
                           max_floor + stiffness.shape[2])) * numpy.nan
    new_mask = numpy.ones((stiffness.shape[0], stiffness.shape[1], 
                           max_floor + stiffness.shape[2]))
    if type(stiffness) == numpy.ma.core.MaskedArray:
        data = stiffness.data
        mask = stiffness.mask
    else:
        data = stiffness
        mask = numpy.zeros_like(stiffness)
        
    max_depth = stiffness.shape[2]
    for x in range(stiffness.shape[0]):
        for y in range(stiffness.shape[1]):
            try:
                last_data = numpy.nonzero(mask[x, y])[0][0]
            except IndexError:
                last_data = len(mask[x, y]) - 1
            begin_pos = floor[x,y]
            end_pos = last_data + begin_pos + 1
            try:
                new_array[x, y, begin_pos:end_pos] = data[x, y][last_data::-1]
            except ValueError as error:
                print "Error in post_proc line 437:"
                print x, y
                print begin_pos, end_pos
                print new_array[x, y, begin_pos:end_pos].shape
                print data[x, y][last_data::-1].shape
                print error
            new_mask[x, y, begin_pos:end_pos] = mask[x, y][last_data::-1]
            
    # Keep only the usefull informations. Array was build too large, but we want
    # to return a well sized array
    first_empty = 0
    nbr_empty = 1
    for depth in range(new_mask.shape[2]):
        if new_mask[:, :, depth].all():
            if depth == first_empty + nbr_empty:
                nbr_empty += 1
            else:
                first_empty = depth
                nbr_empty = 1
    if first_empty == 0 and nbr_empty == 1:
        # This mean there is no empty slice...
        first_empty = new_array.shape[2]
    ##
    ## Stores the result in a matrix that supports the mask, to get rid
    ## of the nan values
    ##
    return numpy.ma.array(new_array[:, :, :first_empty],
                          mask=new_mask[:, :, :first_empty])

def generate_path(pos_a, pos_b):
    """
        Generate a linear path between the point a and the point b.
        
        >>> generate_path([0, 0], [1, 1])
        array([[ 0.,  1.],
               [ 0.,  1.]])
        
        >>> generate_path([0, 0], [1, 5])
        array([[ 0.,  0.,  0.,  1.,  1.,  1.],
               [ 0.,  1.,  2.,  3.,  4.,  5.]])
        
        >>> generate_path([0, 0], [0, 5])
        array([[ 0.,  0.,  0.,  0.,  0.,  0.],
               [ 0.,  1.,  2.,  3.,  4.,  5.]])

        >>> generate_path([0, 0], [5, 0])
        array([[ 0.,  1.,  2.,  3.,  4.,  5.],
               [ 0.,  0.,  0.,  0.,  0.,  0.]])

        >>> generate_path([2, 2], [8, 5])
        array([[ 2.,  3.,  4.,  5.,  6.,  7.,  8.],
               [ 2.,  3.,  3.,  4.,  4.,  5.,  5.]])
    """
    pos_a[0], pos_a[1], pos_b[0], pos_b[1] = \
                int(pos_a[0]), int(pos_a[1]), int(pos_b[0]), int(pos_b[1])
    if pos_b[0] == pos_a[0]:
        if pos_a[1] > pos_b[1]:
            pos_a, pos_b = pos_b, pos_a
        _length = pos_b[1] - pos_a[1] + 1
        _y = range(pos_a[1], pos_b[1] + 1)
        _x = [pos_b[0]] * _length
    elif pos_b[1] == pos_a[1]:
        if pos_a[0] > pos_b[0]:
            pos_a, pos_b = pos_b, pos_a
        _length = pos_b[0] - pos_a[0] + 1
        _x = range(pos_a[0], pos_b[0] + 1)
        _y = [pos_b[1]] * _length
    else:
        _a = float(pos_b[1] - pos_a[1]) / (pos_b[0] - pos_a[0])
        _b = pos_a[1] - _a * pos_a[0]
        
        if abs(pos_b[0] - pos_a[0]) > abs(pos_b[1] - pos_a[1]):
            if pos_a[0] > pos_b[0]:
                pos_a, pos_b = pos_b, pos_a
            # x0 -> x1 is bigger to y0 -> y1
            _x = range(pos_a[0], pos_b[0] + 1)
            _y = [round(_a * item + _b) for item in _x]
        else:
            if pos_a[1] > pos_b[1]:
                pos_a, pos_b = pos_b, pos_a
            # y0 -> y1 is bigger than x0 -> x1
            _y = range(pos_a[1], pos_b[1] + 1)
            _x = [round((item - _b) / _a) for item in _y]
    return numpy.asarray([_x, _y], dtype=float)

def generate_slice(path, array):
    """
        Generate the slice of the array along the path.
        
        >>> array = [range(10)]*10
        >>> path = generate_path([0, 0], [0, 4])
        >>> generate_slice(path, array)
        [0, 1, 2, 3, 4]
        
        >>> path = generate_path([0, 1], [4, 1])
        >>> generate_slice(path, array)
        [1, 1, 1, 1, 1]
    """
    return [array[int(_x)][int(_y)] for _x, _y in zip(path[0], path[1])]

if __name__ == "__main__":
    
    try:
        matplotlib.__version__
    except NameError:
        import matplotlib
        matplotlib.use('Agg')
        from matplotlib import pylab
    import doctest
    doctest.testmod()
