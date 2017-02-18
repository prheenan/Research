#! /usr/bin/env python
#-*- coding: iso-8859-1 -*-
'''
The curve module contains function that are usefull for curve data processing
'''

## Fonctions pour Open Fovea.
#import pdb
import csv
#import warnings
#warnings.simplefilter("ignore")
#warnings.filterwarnings('error', category=Warning)

import numpy as num
#from scipy.optimize import leastsq
from scipy import odr
from copy import deepcopy
CURVE_TEST = '../../docs/data/bugged_curve/curve_bug_1.csv'
CURVE_EVENT = ['../../docs/data/curves/event_1.csv',
               '../../docs/data/curves/event_2.csv',
               '../../docs/data/curves/event_3.csv',
               '../../docs/data/curves/WLC.csv']
FIGURE_FOLDER = 'tests/'


def open_csv_curve(file_name):
    """
        Open the curve in csv file.
    """
    fid = open(file_name)
    temp_file = csv.reader(fid)
    # The two first lines are descripion.
    temp_file.next()
    temp_file.next()
    # Get the datas
    arr_list = [item for item in temp_file]
    curves = [num.asarray(
                    [item[place] for item in arr_list if len(item[place])],
                    dtype=num.float)
             for place in range(len(arr_list[0]))]
    return curves


def find_poc(curve_x, curve_y, threshold=1, method='curve_fit',
             limit_slide=False, len_slice=0.3):
    """
    Finds the point of contact of a force curve.

    >>> curve = open_csv_curve(CURVE_TEST)
    >>> result = find_poc(curve[0], curve[1], threshold=1, method='curve_fit')
    >>> result.keys()
    ['deriv', 'Poly Fit', 'Slice', 'PoC', 'Error']

    curve_x and curve_y = array

    threshold : integer
                how much the curve has to diverge from the fit to be a point
                of contact. (threshold * noise)

    method : string
             Defines the method used to detect the point of contact.

            * 'deriv' (default) : uses the derivative method.
                                When the derivative of the curve is zero (then
                                it's linear), the point of contact is detected.
                                The advantage of this method is to get rid of
                                irregularities observed in the off-contact part
                                of the curve.

            * 'curve_fit' : uses the curve fit method.
                          A linear fit is computed from the off-contact part of
                          the curve. The point of contact is then detected at
                          the point when the curve deviates from this linear
                          fit. The advantage of this method is that it does not
                          depend on the shape of the on contact part of the
                          curve.
    """
    if len(curve_x) <= 2:
        return None
    if curve_x == None or curve_y == None:
        return None
    if method == 'deriv':
        return find_poc_deriv(curve_x, curve_y, threshold, limit_slide,
                              len_slice)
    elif method == 'curve_fit':
        try:
            result = find_poc_curve_fit(curve_x, curve_y, threshold,
                                    limit_slide, len_slice)
        except TypeError as detail:
            if 'expected non-empty vector' in detail.message:
                # This error occurs when a too small curve is given.
                return None
            else:
                raise detail
        else:
            result['deriv'] = None
        return result


def platgliti(curve):
    """
        Return the slice of the first non flat part. This is usefull when a
        curve is complemented with the identical values at the end.
    """

    indice = 0
    prev_val = curve[indice]
    indice += 1
    next_val = curve[indice]
    while prev_val == next_val and indice + 1 < len(curve):
        indice += 1
        next_val, prev_val = curve[indice], next_val
    return indice


def true_data_indice(curve_x):
    """
        Return the true length curve. It removes all similar values that are at
        the beginning or at the end of the curve.
    """
    new_curve_x = curve_x[:-1] - curve_x[1:]
    non_zero = new_curve_x.nonzero()[0]
    start_indice = non_zero[0]
    stop_indice = non_zero[-1] + 2

    return start_indice, stop_indice


def find_poc_curve_fit(curve_x, curve_y, threshold, limit_slide=False,
                       len_slice=0.3):
    '''
    Finds the point of contact, using the method of curve fit.

    This methods makes a first order fits the off contact part of the curve.
    The point of contact is then detected when the force curve deviates from
    the linear fit. The deviation is related to the error measured during the
    fit.
    '''
    loop = {'slope_negative': 1,
            'step': 0}
    reversed_curve = 0
    if curve_y[0] > curve_y[-1]:
        reversed_curve = 1
        curve_y = curve_y[::-1]
        curve_x = curve_x[::-1]
    # Linear fit of the segment slice
    _ind = platgliti(curve_y)
    _true_indice = true_data_indice(curve_x)
    #_true_indice = (0, len(curve_x))
    #print _true_indice
    #print _ind
    _length = _true_indice[1] - _true_indice[0]

    seg_slice = slice(_ind, int(_length * len_slice + _ind))
    # the fit is :
    # p(x) = coeff[0]*x + coeff[1]
    coeff = num.polyfit(curve_x[seg_slice],
                        curve_y[seg_slice], 1, full=True)[0]
    #
    # We make a loop by testing the slope of the fit.
    # If this slope is negative, we look a little bit forward in the curve
    # until the slope becomes positive
    #
    fit = num.polyfit(curve_x[seg_slice], curve_y[seg_slice], 1, full=True)
    # Look in the curve slope was positive at the beginning
    # i.e. :
    # \
    #  \     _.-
    #   \̣.-`
    if fit[0][0] >= 0:
        fall_off = 1
    else:
        fall_off = 0
    while loop['slope_negative'] and (seg_slice.stop + 1 < len(curve_x)):
        new_seg_slice = slice(seg_slice.start + 2, seg_slice.stop + 2)
        fit = num.polyfit(curve_x[new_seg_slice], curve_y[new_seg_slice],
                          1, full=True)
        if coeff[0] < fit[0][0] and fit[0][0] <= 0:
            loop['step'] = 0
            if limit_slide and (seg_slice.stop + 1 >= len(curve_x) / 2):
                # The slice is going farer than the half of the curve. Then,
                # if we choosed to limit the sliding, stop at this step.
                loop['slope_negative'] = 0
            if seg_slice.stop > len(curve_x):
                # The slice is going outside the curve, then stop anyway.
                loop['slope_negative'] = 0
        elif fall_off and fit[0][0] >= 0:
            loop['step'] = 0
        elif fall_off and fit[0][0] < 0:
            # we directly stop here and go back
            loop['slope_negative'] = 0
            [new_seg_slice, fit] = go_back(curve_x, curve_y, seg_slice)
#            new_seg_slice = slice(seg_slice.start - 4, seg_slice.stop - 4)
#            fit = num.polyfit(curve_x[new_seg_slice], curve_y[new_seg_slice],
#                              1, full = True)
        else:
            loop['step'] = loop['step'] + 1
        coeff = fit[0]
        if loop['step'] > 2:
            coeff[0] = fit[0][0]
            coeff[1] = fit[0][1]
            loop['slope_negative'] = 0
        seg_slice = new_seg_slice
    #
    # We take the selected slice to detect the point of contact
    #
    error = (num.sqrt(fit[1]) / 5) * threshold
    # Find when the curve goes away from the fit
    #pdb.set_trace()
    if not len(error):
        error = num.array([0])
    _ind = num.polyval([coeff[0], coeff[1] + error], curve_x)
    point_of_contact = 0

    for item_nbr in xrange(len(curve_y)):
        if _ind[item_nbr] >= curve_y[item_nbr]:
            point_of_contact = item_nbr

    if reversed_curve:
        curve_y = curve_y[::-1]
        point_of_contact = len(curve_y) - point_of_contact - 1
        seg_slice = slice(len(curve_y) - seg_slice.stop,
                          len(curve_y) - seg_slice.start)
    return {'PoC': point_of_contact,
            'Poly Fit': [coeff[0], coeff[1]],
            'Error': error,
            'Slice': seg_slice}


def go_back(curve_x, curve_y, seg_slice):
    """
        Go back in the fit to find the best place to find poc.
    """
    fit = num.polyfit(curve_x[seg_slice], curve_y[seg_slice], 1,
                      full=True)
    better = 0
    while not better:
        new_seg_slice = slice(seg_slice.start - 2, seg_slice.stop - 2)
        new_fit = num.polyfit(curve_x[seg_slice], curve_y[seg_slice], 1,
                      full=True)
        if new_fit[1] > fit[1]:
            better = 1
        else:
            seg_slice = new_seg_slice
            fit = new_fit
        if seg_slice.start <= 2:
            better = 1
    return seg_slice, fit


def compute_derivative(curve_x, curve_y, delta=10):
    '''
        Computes the derivative of a curve
    '''
    _my = num.zeros_like(curve_y)
    _mx = num.zeros_like(curve_y)
    for item in range(len(curve_y)):
        t_slice = slice(item, item + delta)
        _my[item] = curve_y[t_slice].mean()
        _mx[item] = curve_x[t_slice].mean()
    d_y = (_my[1:] - _my[:-1])
    return[curve_x[:-1], d_y]


def find_poc_deriv(curve_x, curve_y, threshold=1, limit_slide=False,
    len_slice=0.3):
    '''
        Finds the point of contact, using the derivative of the curve.
    '''
    seg = {'slice': None,  # Contains the slice of the segment
            'coef': None}  # Contains the coefficient of the fit
    # We slips on the curve. Then, we define the same for the next portion of
    # curve
    next_seg = {'slice': None,  # Contains the slice of the segment
            'coef': None}  # Contains the coefficient of the fit

    deriv = compute_derivative(curve_x, curve_y, delta=10)
    # define the slice to compute the fit
    seg['slice'] = slice(int(len(deriv[1]) - len(deriv[1]) * len_slice) - 1,
                      len(deriv[1]) - 1)
    # and compute the fit...
    seg['coef'] = num.polyfit(deriv[0][seg['slice']], deriv[1][seg['slice']],
                           1, full=True)
    # make the slice slips on the curve, untill we find the one with the lower
    # error.
    # Indeed, some curves begin with a big noise that prevent a good point of
    # contact detection.
    loop = {'ok': 0,
            'i': 0,
            'step': 0}
    while not loop['ok']:
        next_seg['slice'] = slice(seg['slice'].start - loop['step'],
                               seg['slice'].stop - loop['step'])
        next_seg['coef'] = num.polyfit(deriv[0][next_seg['slice']],
                                    deriv[1][next_seg['slice']],
                                    1, full=True)
        if seg['coef'][1] <= next_seg['coef'][1]:
            loop['ok'] = 1
        #else:
        #    i = 0
        if limit_slide and (next_seg['slice'].start < len(deriv[1]) / 2):
            loop['ok'] = 1
        if next_seg['slice'].start < loop['step']:
            loop['ok'] = 1
        seg['coef'] = next_seg['coef']
        loop['i'] += 1
        if loop['i'] > 10:
            seg['slice'] = next_seg['slice']
            loop['ok'] = 1
    # Now, we apply the threshold on the error and find the first point that
    # deviates from the fit (same idea as the old curve_fit algorythm)
    error = seg['coef'][1] * threshold
    seg['coef'] = seg['coef'][0]
    polynom = num.polyval(seg['coef'], deriv[0]) - error
    on_contact = (deriv[1] <= polynom) - 1  # True are the points that are
                                            # close to the fit
    indice = num.asarray(num.nonzero(on_contact)[0])
    # Sometimes, in the on-contact part of the curve, the tips slips on the
    # surface that is indented. The resulting curve looks like this :
    # \
    #  \  /\
    #   `v  \
    #        "_
    #          `¬___________________
    #    |      |
    #    V      V
    # 000100000011111111111111111111 Resulting on_contact
    # 123456789012345678901234567890
    # [4, 11, 12, 13, 14, 15, 16, 17, ...] Resulting indice
    # [7, 1, 1, 1, 1, 1, 1, ...] Resulting diff_indice
    #
    # With the diff_indice, we can then detect such holes. I arbitrary choose
    # 40. It could be nice to set this as a parameter. Will see later.
    diff_indice = indice[1:] - indice[:-1]
    hole = num.nonzero(diff_indice > 40)[0]
    if len(hole):
        i_poc = hole[-0] + 1
    else:
        i_poc = 0
    non_zero_list = num.nonzero(on_contact)[0]
    if len(non_zero_list):
        point_of_contact = num.nonzero(on_contact)[0][i_poc]
    else:
        point_of_contact = 0
    return {'PoC': point_of_contact,
            'Poly Fit': [0, seg['coef'][0]],
            'Error': error,
            'Slice': seg['slice'],
            'deriv': deriv}


def compute_indentation(curve_x, curve_y,
                        deflection_sensitivity, spring_constant):
    '''
    Compute the indentation curve.

    * Parameters :
        curve_x : 1d array
            Values of the x coordinate.
        curve_y : 1d array
            Values of the y coordinates.
        deflection_sensitivity : float
            Parameter of the glass indentation.
        spring_constant : float
            sping constant of the cantilever.

    * Returns :
        indentation : 1d array
            Value of the x coordinate.
        force : 1d array
            Value of the y coordinate.
    '''
    #deflection_sensitivity=1
    if curve_x == None or curve_y == None or len(curve_x) <= 1:
        return [None, None]
    indice = true_data_indice(curve_x)
    length = len(curve_x)

    curve_x = curve_x[indice[0]:indice[1]]
    curve_y = curve_y[indice[0]:indice[1]]
    if curve_x[0] < curve_x[-1]:
        indent_curve = curve_x[::-1]
    else:
        indent_curve = curve_x

    indent_curve = indent_curve + num.polyval([deflection_sensitivity, 0],
                                              curve_y)
    force_curve = curve_y * spring_constant

    force_curve = num.r_[force_curve,
                         [force_curve[-1]] * (length - len(force_curve))]
    indent_curve = num.r_[indent_curve,
                          [indent_curve[-1]] * (length - len(indent_curve))]
    return [indent_curve, force_curve]


def correct_curve_drop(curve_y):
    '''
    Corrects the drops that occurs at the beginning of some force-distance
    curves. This function corrects a bug from Veeco AFM.
    '''
    first_non_zero = curve_y.nonzero()[0][1] + 1
    for indice in range(first_non_zero, len(curve_y)):
        curve_y[indice] = curve_y[first_non_zero - 1]


def segment_curve(curve_x, curve_y, segment_number, segment_deep):
    """
    Segmentation of curve.

    Parameters : 
        * curve_x : 1d array
                Values of the x coordinate.
        * curve_y : 1d array
                Values of the y coordinate.
        * segment_number : int
                specifies the number of segments to generate
                -1 : create one segment from the begining to
                     the end of the curve.
                0 : create as many segment as possible.
                >0 : specifies the exact number of segment. If the
                     curve is smaller, less segment is created.
        * segment_deep : int
                specifies the deep (curve_x value) of each segment.
                in case segment_number = -1, this parameter is not
                considered.

    Returns :
        * parts_x : 1d array
              segment_number+1 array containing the curve_x values that
              correspond to each segment
        * parts_y : 1d array
              segment_number+1 array containing the curve_y values that
              correspond to each segment

        for example : segment n is made of parts_x[n-1],parts_x[n]
                      and parts_y[n-1],parts_y[n]
    """
    ## Initialisation des valeurs
    ##
    #
    if not len(curve_y):
        if segment_number == 0:
            return [[num.nan], [num.nan]]
        else:
            return [num.zeros(segment_number) * num.nan,
                    num.zeros(segment_number) * num.nan]
    if curve_y[0] > curve_y[-1]:
        curve_x = curve_x[::-1]  # We invert the curve to be more convenient
        curve_y = curve_y[::-1]

    meta_parts_x = []
    meta_parts_y = []

    if segment_number == -1:
        # We create one segment that correspond to the whole curve.
        meta_parts_y.append(curve_y)
        meta_parts_x.append(curve_x)
        return [meta_parts_x, meta_parts_y]
    elif segment_number == 0:
        segment_number = num.inf
    if len(curve_x):
        curve_size = curve_x.max() - curve_x.min()
    else:
        curve_size = 0
    size_done = segment_deep
    nbr_selection = 1
    ## Element zero
    parts_x = []
    parts_y = []
    if curve_size:
        parts_x.append(curve_x.min())
        parts_y.append(curve_y[num.nonzero(curve_x == parts_x[0])[0][0]])
    else:
        parts_x.append(num.nan)
        parts_y.append(num.nan)
    ## From the first to the last segment
    all_indexes = [0]
    while curve_size > size_done and nbr_selection <= segment_number:
        parts_x.append(parts_x[-1] + segment_deep)
        # Let's find the next point
        indexes = num.nonzero(curve_x >= parts_x[nbr_selection])[0]
        indexes = indexes[0]
        if indexes == 0:
            indexes = 1
        all_indexes.append(indexes)
        # The point stands between index-1 and index. We then make a linear fit
        # and we'll find the value. Done between index-1 and index+1 'cause it
        # does _not_ take the last value.
        polynome = num.polyfit(curve_x[indexes - 1:indexes + 1],
                               curve_y[indexes - 1:indexes + 1],
                               1)
        parts_y.append(num.polyval(polynome, parts_x[nbr_selection]))
        this_slice = slice(all_indexes[nbr_selection - 1],
                           all_indexes[nbr_selection])
        if parts_x[-2] != curve_x[this_slice.start]:
            tmp_x = [parts_x[-2]]
            tmp_y = [parts_y[-2]]
        else:
            tmp_x = []
            tmp_y = []
        tmp_x.extend(curve_x[this_slice])
        tmp_x.append(parts_x[-1])
        #print tmp_x
        #tmp_y = [parts_x[-1]]
        tmp_y.extend(curve_y[this_slice])
        tmp_y.append(parts_y[-1])
        meta_parts_x.append(tmp_x)
        meta_parts_y.append(tmp_y)
        nbr_selection = nbr_selection + 1
        size_done = size_done + segment_deep
    if nbr_selection <= segment_number and segment_number is not num.inf:
        ## adding NAN values
        while nbr_selection <= segment_number:
            meta_parts_x.append([num.nan])
            meta_parts_y.append([num.nan])
            nbr_selection += 1
    return [meta_parts_x, meta_parts_y]


def compute_stiffness(indent_curve, force_curve, model='Sphere', tip_carac=40,
                      poisson_ratio=0.3, method='Raw'):
    """
    Compute : the stiffness of the different depth.

    indent_curve and force_curve : the verctors return by segment_curve.

    model : describe the Hertz model used and can be 'Sphere' or 'Cone'.

    tip_carac : the characteristic of the tip. Radius and Semi-opening angle in
            radian for Sphere or Cone model respectively.
    poisson_ratio : the poisson ratio of the indented substrate.
                    Common values : 0.3 for normal cells
                                    0.5 for rubber
    method : string
             'Raw' for a fit of the raw data with the model.
             'Linear' for a linearized fit.
             'Extrema' for a linearized fit from the extrem points.
             'Median' for a computation of the ym from the median point in the
                      segment.
    """
    if indent_curve == None or force_curve == None:
        return num.ones(1) * num.nan
    ## The indentation and force have to begin from 0...
    stop = False
    if not len(indent_curve):
        stop = True
    elif type(indent_curve[0]) == list:
        stop = len(indent_curve[0]) == 1 and num.isnan(indent_curve[0][0])
    elif type(indent_curve[0]) in [num.float64, float]:
        stop = num.isnan(indent_curve[0])
    if stop:
        return num.ones(len(indent_curve)) * num.nan
    # Make the point of contact being at zero:
    zero_point = indent_curve[0][0]
    indent_curve = [item - zero_point for item in indent_curve]
    zero_point = force_curve[0][0]
    force_curve = [item - zero_point for item in force_curve]
    # Now, we prepare the fit
    young_modulus = []
    for indent, force in zip(indent_curve, force_curve):
        if num.isnan(indent[0]):
            young_modulus.append(num.nan)
        elif method == 'Raw':
            # Choose the initial parameters :
            p = [100]

            def ym_model(p, indent):
                if model == 'Sphere':
                    force = 4 / 3. * p[0] / (1 - poisson_ratio ** 2) *\
                            num.sqrt(tip_carac) * indent ** 1.5
                elif model == 'Cone':
                    force = 2. / num.pi * p[0] / (1 - poisson_ratio ** 2) *\
                            num.tan(tip_carac) * indent ** 2
                return force

            try:
                fit_it = odr.ODR(odr.Data(indent, force),
                                 odr.Model(ym_model), p)
                # And we run it
                out = fit_it.run()
                # Get the result and rescale it.
                young_modulus.append(out.beta[0])
                if out.beta[0] == 100:
                    print indent
                    out.pprint()
            except IndexError, message:
                if type(indent) in [float, num.float64] and num.isnan(indent):
                    young_modulus.append(num.nan)
                else:
                    raise IndexError(message)
        elif method == 'Linear':
            # Choose the initial parameters :
            p = [100]

            def ym_model(p, d_indent):
                #slope = (force[1:] - force[:-1]) /(indent[1:] - indent[:-1])
                if model == 'Sphere':
                    d_force = 4 / 3. * p[0] / (1 - poisson_ratio ** 2) *\
                            num.sqrt(tip_carac) * d_indent
                elif model == 'Cone':
                    d_force = 2. / num.pi * p[0] / (1 - poisson_ratio ** 2) *\
                            num.tan(tip_carac) * d_indent
                return d_force
            try:
                if model == 'Sphere':
                    indent = indent ** (3 / 2.)
                elif model == 'Cone':
                    indent = indent ** 2
                d_force = num.diff(force)
                d_indent = num.diff(indent)
                fit_it = odr.ODR(odr.Data(d_indent, d_force),
                                 odr.Model(ym_model), p)
                # And we run it
                out = fit_it.run()
                # Get the result and rescale it.
                young_modulus.append(out.beta[0])
            except IndexError, message:
                if type(indent) in [float, num.float64] and num.isnan(indent):
                    young_modulus.append(num.nan)
                else:
                    raise IndexError(message)
        elif method == 'Extrema':
            # Modification of the axes ...
            if model == 'Sphere':
                indent = indent ** (3 / 2.)
            elif model == 'Cone':
                indent = indent ** 2
            # Compute slopes
            try:
                slope = (force[-1] - force[0]) / (indent[-1] - indent[0])
            except IndexError, message:
                if type(indent) in [float, num.float64] and num.isnan(indent):
                    young_modulus.append(num.nan)
                else:
                    raise IndexError(message)
            finally:
                ## Compute the young modulus according to the models
                if model == 'Sphere':
                    young_modulus.append((1 - poisson_ratio ** 2) *
                                         slope *
                                         3 / (4 * num.sqrt(tip_carac)))
                elif model == 'Cone':
                    young_modulus.append((1 - poisson_ratio ** 2) *
                                         slope *
                                         num.pi / (2 * num.tan(tip_carac)))
        elif method == 'Median':
            # Modification of the axes ...
            if model == 'Sphere':
                indent = indent ** (3 / 2.)
            elif model == 'Cone':
                indent = indent ** 2
            m_force = num.median(force)
            m_indent = num.median(indent)
            if model == 'Sphere':
                young_modulus.append((1 - poisson_ratio ** 2) *
                                         (m_force / m_indent) *
                                         3 / (4 * num.sqrt(tip_carac)))
            elif model == 'Cone':
                young_modulus.append((1 - poisson_ratio ** 2) *
                                         (m_force / m_indent) *
                                         num.pi / (2 * num.tan(tip_carac)))
    return num.array(young_modulus)


def event_find(curve_x, curve_y, deflection_sensitivity, spring_constant,
               weight=2.5, fit_model=None, poc=0, baseline=None, debug=0):
    '''
        This function finds the protein - protein unbinding events from a
        retraction curve.

        * Parameters :
            curve_x : 1d array
                Values of the x coordinate.
            curve_y : 1d array
                Values of the y coordinates.
            deflection_sensitivity : float
                The deflection sensitivity is used to compute properties of the
                event.
            spring_constant : float
                The property of the cantilever to compute the force of the
                event.
            weight : float
                Threshold to detect the event. Smallest is more sensitive.
            fit_model : str
                'wlc' for worm like chain model
                'fjc' for freely jointed chain model
            poc : int
                The index of the point of contact in the curve. Used for the
                fit.
            baseline : [a, b] -> y = ax + b
                The indices of the baseline of the off contact part of the
                curve. Used to compute the fit.

        * Returns :
            event_list : list
                List of events as returned by event_properties.
    '''
    if curve_x == None or curve_y == None:
        return None
    # find curve jumps...
    [curve_pos, curve_mean, curve_std] = find_curve_jump(curve_y)
    indent = compute_indentation(curve_x, curve_y, deflection_sensitivity,
                                   spring_constant)
    distances = -indent[0] + max(indent[0])
    threshold = weight * curve_pos * curve_std
    positions = threshold > 1
    event_list = None
    # Find the point of contac
    poc_xy = [curve_x[poc], curve_y[poc]]
    if baseline is not None:
        f_y = baseline[0] * curve_x + baseline[1]
        try:
            new_poc = num.nonzero(curve_y < f_y)[0][0] - 1
        except:
            new_poc = poc
        poc_xy = [curve_x[new_poc], curve_y[new_poc]]
    if sum(positions):
        event_list = event_complete(curve_x, curve_y, positions)
        for event in event_list:
            event = event_properties(curve_x, curve_y, event, spring_constant,
                                 distances, fit_model, poc_xy=poc_xy)
    if not debug:
        return event_list
    else:
        return [event_list, threshold]


def event_complete(curve_x, curve_y, position):
    '''
    from a vector which points the putative events position, this function
    extracts and sort and adds some informations like force,...
    '''
    point_list = position.nonzero()[0]
    # Finds adjascent points that defines in fact a single event...
    event_list = list()
    first_pos = 0
    to_remove = list()
    for i in xrange(len(point_list) - 1):
        if point_list[i] + 1 != point_list[i + 1]:
            event_list.append({'X': point_list[i],
                               'Slice': slice(point_list[first_pos],
                                              point_list[i] + 1)})
            first_pos = i + 1
        if i + 1 == len(point_list) - 1:
            event_list.append({'X': point_list[i],
                               'Slice': slice(point_list[first_pos],
                                              point_list[i + 1] + 1)})
            first_pos = i + 1
    if len(point_list) == 1:
        event_list = [{'X': point_list[0],
                       'Slice': slice(point_list[0],
                                      point_list[0] + 1)}]
    # Complete all the events.
    for event_nbr in xrange(len(event_list)):
        event = event_list[event_nbr]
        # Complete the events to the right
        next_is_higher = 1
        right_pos = event['Slice'].stop
        while next_is_higher and right_pos + 1 < len(curve_y):
            if curve_y[right_pos + 1] >= curve_y[right_pos]:
                right_pos = right_pos + 1
            else:
                right_pos = right_pos + 1
                next_is_higher -= 1
        # And complete to the left.
        previous_is_lower = 1
        left_pos = event['Slice'].start
        while previous_is_lower:
            if curve_y[left_pos - 1] < curve_y[left_pos]:
                left_pos = left_pos - 1
                #previous_is_lower = 1
            elif curve_y[left_pos - 1] < curve_y[left_pos + 1]:
                left_pos = left_pos - 1
                #previous_is_lower = 1
            else:
                #print curve_y[left_pos]
                #print curve_y[left_pos + 1]
                left_pos = left_pos - 1
                previous_is_lower -= 1
        # We are at the bottom of the event. Let's see the event jump...
        nbr_test = round(len(curve_y) / 130)
        previous_is_higher = nbr_test
        while previous_is_higher and left_pos > 0:
            if curve_y[left_pos - 1] >= curve_y[left_pos]:
                left_pos = left_pos - 1
#                previous_is_higher = 1
            elif curve_y[left_pos - 1] >= curve_y[left_pos + 3]:
                left_pos = left_pos - 1
                #previous_is_lower =
            else:
                left_pos = left_pos - 1
                previous_is_higher = previous_is_higher - 1
        left_pos = left_pos + nbr_test
        #print curve_x[left_pos]
        if left_pos <= nbr_test or curve_x[left_pos] <= 5:
            # This is a jump off contact
            to_remove.append(event_nbr)
        event['Slice'] = slice(left_pos, right_pos + 1)
    for rm_nb in range(len(to_remove)):
        #print(to_remove[rm_nb] - rm_nb)
        event_list.pop(to_remove[rm_nb] - rm_nb)
    return event_list


def event_properties(curve_x, curve_y, event, spring_constant, distances,
                     fit_method=None, poc_xy=None):
    '''
    Given an event, event_properties computes the properties of the events. The
    result is returned as a dictionnary with the following keys :
    'min' : the minimum of the event, just before the bound break
    'max' : the maximum of the event, just after the bound break
    'slope' : the slope of the curve that defins the bound breaking
    'force' : the unbinding force of the event
    'loading_rate' : the loading rate of the event
    'fit_length' : the length of the protein found via the fit.
    'fit_plength' : the persistent length of the protein.
    '''
    # Finds minimum
    value = min(curve_y[event['Slice']])
    event['min'] = int(list(curve_y[event['Slice']]).index(value) + \
                   event['Slice'].start)
    if event['min'] - event['Slice'].start == 1:
        # Minimum and start are too close. Increase the slice size of the event
        # to the left.
        event['Slice'] = slice(event['Slice'].start - 1, event['Slice'].stop)
    # Finds maximum
    value = max(curve_y[event['min']:event['Slice'].stop])
    event['max'] = list(curve_y[event['min']:
                        event['Slice'].stop]).index(value) + event['min']
    if event['max'] > len(curve_x):
        event['max'] = len(curve_x)
    # Find jumpSlope
    event['slope'] = ((curve_y[event['max']] - curve_y[event['min']]) /
                       (curve_x[event['max']] - curve_x[event['min']]))
    # Find loading rage
    if event['Slice'].start == event['min']:
        event['loading_rate'] = num.nan
    else:
        event['loading_rate'] = -((curve_y[event['Slice'].start] -
                                                    curve_y[event['min']]) /
                              (curve_x[event['Slice'].start] -
                                                    curve_x[event['min']]))
    # Find force
    event['force'] = (curve_y[event['max']] - curve_y[event['min']]) * \
                     spring_constant
    event['force_base'] = (curve_y[-1] - curve_y[event['min']]) * \
                          spring_constant
    #print('Force : %f, Force_base : %f')%(event['force'], event['force_base'])
    #print('Curve 0 : %f, Curve -1 : %f, ')%(curve_y[0], curve_y[-1])
    #print('Event max : %f')%(curve_y[event['max']])
    # find distance
    event['dist'] = distances[event['min']]
    ## integrate with PoCposition
    # Perfor the fit on the event.
    if event_fit:
        [length, plength, fit_x, fit_y] = event_fit(curve_x, curve_y, event,
                                            method=fit_method,
                                            spring_constant=spring_constant,
                                            poc_xy=poc_xy)
    else:
        [length, plength, fit_x, fit_y] = [None, None, None, None]

    event['fit_length'] = length
    event['fit_plength'] = plength
    event['fit_x'] = fit_x
    event['fit_y'] = fit_y
    return event


def find_curve_jump(curve):
    '''
    This function finds jumps in a curve. It returns three vectors of the same
    size of the input vector.
    curve_jump : filled with zeros and ones. Ones depicted position of jumps.
    curve_mean : contain the flatten curve
    curve_std  : contain the standard deviation measured along the curve.
    '''
    window_size = 10
    curve_jump = num.zeros(len(curve))
    curve_mean = deepcopy(curve)
    curve_std = num.zeros(len(curve))
    portion_std = num.std(curve)
    for point_nbr in xrange(len(curve) - window_size):
        curve_portion = curve[point_nbr:point_nbr + window_size]
        indice = 0
        middle_point = curve_portion[indice]
        portion_mean = num.mean(curve_portion)
        portion_std = num.std(curve_portion)
        curve_mean[point_nbr] = portion_mean
        curve_std[point_nbr] = portion_std
        if middle_point < portion_mean - portion_std * 2:
            if point_nbr + indice:
                curve_jump[point_nbr + indice] = \
                curve_jump[point_nbr + indice - 1] + 1
            else:
                curve_jump[point_nbr + indice] = 1
    return curve_jump, curve_mean, curve_std


def event_fit(curve_x, curve_y, event, method,
              spring_constant=0.06, poc_xy=None):
    '''
        Fit the event with the correct method.

        * Parameters :
            curve_x : 1d array
                Values of the x coordinate.
            curve_y : 1d array
                Values of the y coordinates.
            event : dict
                dictionnary as returned by event_find function.
            method : str
                'wlc' for worm like chain model
                'fjc' for freely jointed chain model
            spring_constant : float
                The spring constant of the used cantilever
            poc : int
                The indice of the point of contact as detected in the approach
                curve.

        * Returns :
            length : 1d vector
                length[0] = length in nm.
                lenght[1] = error in nm.
            plength : 1d vector
                plength[0] = persistent length in nm.
                plenght[1] = error in nm.
            fit_x : 1d array
                contains the x values of the fit on the curve.
            fit_y : 1d array
                contains the y values of the fit on the curve.

        As an example, we can get a curve with event and find event with
        event_find :

        >>> curve = open_csv_curve(CURVE_EVENT[0])
        >>> events = event_find(curve[2], curve[3], 0, 0.06)

        In order to fit with the worm like chain model :

        >>> [length, plength, fit_x, fit_y] = \
                event_fit(curve[2], curve[3], events[1], method='wlc',\
                spring_constant=0.0554)

        The result are the length and the persistent length :

        >>> print('Length : %.5f, error : %.5f')%(length[0], length[1])
        Length : 282.42567, error : 1.19537
        >>> print('Plength : %.5f, error : %.5f')%(plength[0], plength[1])
        Plength : 0.26088, error : 0.01271

        Let's plot the result :

        >>> import pylab
        >>> wlc_plot = pylab.plot(curve[2], curve[3])
        >>> wlc_plot = pylab.plot(fit_x, fit_y)
        >>> pylab.savefig(FIGURE_FOLDER + 'wlc_plot.png')

        Here is the same for the free joint chain model :

        >>> [length, plength, fit_x, fit_y] = \
                event_fit(curve[2], curve[3], events[1], method='fjc',\
                spring_constant=0.0554)

        You can see, in that case, that the result are close to the wlc model :

        >>> print('%.3f, %.3f')%(length[0],length[1])
        263.263, 0.865
        >>> print('%.3f, %.3f')%(plength[0],plength[1])
        0.253, 0.010

        Let's plot the result :

        >>> pylab.cla()
        >>> fjc_plot = pylab.plot(curve[2], curve[3])
        >>> fjc_plot = pylab.plot(fit_x, fit_y)
        >>> pylab.savefig(FIGURE_FOLDER + 'fjc_plot.png')
    '''
    # extract the portion of the curve that contains the event
    ar_event = [curve_x[event['Slice']], curve_y[event['Slice']]]
    if poc_xy is None:
        poc_xy = [curve_x[0], None]
        #max_y = None
    if method == 'wlc':
        result, fit, p_eval, maxy = __wlc_fit(ar_event[0], ar_event[1],
                                    poc=poc_xy[0],
                                    spring_constant=spring_constant,
                                    max_y=poc_xy[1])
        # prepare the array for the fit
        poc_indice = num.nonzero(curve_x >= poc_xy[0])[0][0]
        new_slice = slice(poc_indice, event['Slice'].stop)  # event['min'])
        fit_x = curve_x[new_slice]
        fit_y = p_eval(result.beta, fit_x)
    elif method == 'fjc':
        result, fit, p_eval, maxy = __fjc_fit(ar_event[0], ar_event[1],
                                    poc=poc_xy[0],
                                    spring_constant=spring_constant,
                                    max_y=poc_xy[1])
#        # prepare the array for the fit
        min_y = min(ar_event[1])
#        if baseline is not None:
#            max_y = baseline_y[0]
#        else:
#            max_y = max(ar_event[1])
        fit_y = num.arange(min_y, maxy, (maxy - min_y) / 100)
        fit_x = p_eval(result.beta, fit_y) + poc_xy[0]
        #fit_x, fit_y = None, None
    elif method in [None, 'None']:
        return [None, None, None, None]
    else:
        raise TypeError('Method %s is not known.' % method)
    return [[result.beta[0] * 1e9, result.sd_beta[0] * 1e9],
            [result.beta[1] * 1e9, result.sd_beta[1] * 1e9],
            fit_x, fit_y]


def __wlc_fit(x, y, poc=0, spring_constant=0.06, baseline_y=None, max_y=None):
    '''
        Perform Worm Like Chain fit on the curve portion.

        It's better to use "event_fit" with the option "method=wlc" as it
        returns more comprehensive results.

        * Parameters :
            curve_x : 1d array
                Values of the x coordinate.
            curve_y : 1d array
                Values of the y coordinates.
            poc : int
                The index of the point of contact.
            spring_constant : float
                The spring constant of the cantilever (in N/m)

        * Returns :
            out : odr.ODR object
                The object from which you can get the parameters of the fit.
            result : The parameter of the fit # we can remove it...

            pevalr : fct
                function to evaluate the fit.
            rec_maxy
                A conversion value to evaluate the fit.
    '''
    # First, we take the beginning to the lowest point.
    new_x, new_y, rec_maxy = __prepare_event_for_fit(x, y, spring_constant,
                                                     poc=poc, max_y=max_y)
    # The initial condition
    T = 285
    L = new_x[-1] + (new_x[-1] * 0.1)
    Lp = 0.16e-9
    p0 = [L, Lp]

    def peval(p, x):
        """
            Define the function to fit
        """
        # p[0] = L
        # p[1] = Lp
        Kb = 1.381e-23  # Boltzmann constant

        old_settings = num.seterr(all='ignore')
        _peval_result = ((Kb * T / p[1]) *
                    (((x) / p[0]) + 1 / (4 * (1 - (x) / p[0]) ** 2) - 1 / 4))
        num.seterr(**old_settings)  # restore settings
        return _peval_result

    def pevalr(p, x):
        """
            This is the function that will be return in order to complete the
            fit. Usefull for pretty plots.
        """
        # Lp = 1e-9
        new_x = (x - poc) * 1e-9
        result = peval(p, new_x)
        result = -result * 1e9 / spring_constant + rec_maxy
        return result

    # Now, we prepare the fit
    fit_it = odr.ODR(odr.Data(new_x, new_y), odr.Model(peval), p0)
    # And we run it
    out = fit_it.run()
    # Get the result and rescale it.
    result = [new_x, peval(out.beta, new_x)]
    result[0] = result[0] * 1e9 + poc
    result[1] = -(result[1] * 1e9 / spring_constant) + rec_maxy
    return [out, result, pevalr, rec_maxy]


def __prepare_event_for_fit(x, y, spring_constant, poc=0, max_y=None):
    min_index = num.where(y == min(y))[0][0]
    if max_y is not None:
        rec_maxy = max_y
    else:
        rec_maxy = max(y)
    new_x = x[0:min_index + 1] - poc  # - x[poc]
    new_y = y[0:min_index + 1]
    # Then, we rescale to be in the metric system (we were in nm)
    new_x = (new_x) * 1e-9
    new_y = -(new_y - rec_maxy) * 1e-9 * spring_constant  # * 0.06
    return new_x, new_y, rec_maxy


def __fjc_fit(x, y, poc=0, spring_constant=0.06, max_y=None):
    '''
        Perform Worm Like Chain fit on the curve portion.

        It's better to use "event_fit" with the option "method=fjc" as it
        returns more comprehensive results.

        * Parameters :
            curve_x : 1d array
                Values of the x coordinate.
            curve_y : 1d array
                Values of the y coordinates.
            poc : int
                The index of the point of contact.
            spring_constant : float
                The spring constant of the cantilever (in N/m)

        * Returns :
            out : odr.ODR object
                The object from which you can get the parameters of the fit.
            result : The parameter of the fit # we can remove it...

            pevalr : fct
                function to evaluate the fit.
            rec_maxy
                A conversion value to evaluate the fit.
    '''
    # First, we get only the part of the curve that should be fitted.
    new_x, new_y, rec_maxy = __prepare_event_for_fit(x, y, spring_constant,
                                                     poc, max_y=max_y)
    # The initial condition
    T = 285
    L = new_x[-1] + (new_x[-1] * 0.1)
    Lp = 0.4e-9
    p0 = [L, Lp]

    def peval(p, y):
        # p[0] = L
        # p[1] = Lp
        Kb = 1.381e-23  # Boltzmann constant
        # Avoid zero elements in y. As force should be positive, if zero is
        # given, it's the minimum one. Then replace by the non-zero minimum
        # value, which should be little bit higher value.
        y[num.nonzero(y == 0)] = min(y[num.nonzero(y != 0)])
        return (p[0] * (1 / num.tanh((y * p[1]) / (Kb * T)) -
                            (Kb * T) / (y * p[1])))

    def pevalr(p, y):
        """
            This is the function that will be return in order to complete the
            fit. Usefull for pretty plots.
        """
        new_y = -(y - rec_maxy) * 1e-9 * spring_constant
        result = peval(p, new_y)
        result = result * 1e9
        return result
    # Prepare the fit
    fit_it = odr.ODR(odr.Data(new_y, new_x), odr.Model(peval), p0)
    # And we run it
    out = fit_it.run()
    # Get the result and rescale it.
    result = [new_y, peval(out.beta, new_y)]
    result[0], result[1] = result[1], result[0]
    result[0] = result[0] * 1e9
    result[1] = -(result[1] * 1e9 / 0.06 - rec_maxy)
    return [out, result, pevalr, rec_maxy]


if __name__ == '__main__':
    import doctest
    doctest.testmod()
