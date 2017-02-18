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
import os
import csv

import numpy
import urllib2
import datetime

def remove_from_mask(data, mask):
    """
        Remove data from the mask.
        
        Let's have a data of this kind :
        
        >>> import numpy
        >>> data = [{'pos_x' : 0, 'pos_y' : 0, 'item' : 0},\
                    {'pos_x' : 0, 'pos_y' : 1, 'item' : 1},\
                    {'pos_x' : 0, 'pos_y' : 2, 'item' : 2},\
                    {'pos_x' : 0, 'pos_y' : 3, 'item' : 3},\
                    {'pos_x' : 1, 'pos_y' : 3, 'item' : 4},\
                    {'pos_x' : 2, 'pos_y' : 0, 'item' : 5},\
                    {'pos_x' : 2, 'pos_y' : 3, 'item' : 6},\
                    ]
        
        And let's have a mask of this kind :
        
        >>> mask = numpy.array([[True, True, False, False],\
                    [True, True, False, False],\
                    [True, True, False, False],\
                    [True, True, False, False]])
        
        Now, let's apply the remove from mask :
        
        >>> new_list = remove_from_mask(data, mask)
    """
    # Data can be from event structure.
    # In the data structure, we have info about the position on the array. We
    # then exclude data that are in the mask. data['pos_x'] and data['pos_y']
    new_list = []
    for item in data:
        if not mask[item['pos_x'], item['pos_y']]:
            new_list.append(item)
    return new_list

def generate_mask(array, value):
    mask = array < mdata
    return mask

def mask_bool(mask_1, mask_2, fct):
    """
    Make boolean operation on mask arrays.

    Boolean operation can be AND, OR, XOR, NOT

    >>> array_1 = [True, True, False, False]
    >>> array_2 = [False, True, False, True]
    >>> mask_bool(array_1, array_2, "AND")
    array([False,  True, False, False], dtype=bool)
    >>> mask_bool(array_1, array_2, "OR")
    array([ True,  True, False,  True], dtype=bool)
    >>> mask_bool(array_1, array_2, "XOR")
    array([ True, False, False,  True], dtype=bool)
    """

    # If the two arrays does not have the same shape, raise exception
    #if numpy.shape(mask_1) != numpy.shape(mask_2):
    #    raise ValueError, "operands could not be broadcast together with shapes {0} {1}".format(numpy.shape(mask_1), numpy.shape(mask_2))
    #__mask_1 = check_bool(mask_1)
    #__mask_2 = check_bool(mask_2)
    if fct == "AND":
        return numpy.logical_and(mask_1, mask_2)
    elif fct == "OR":
        return numpy.logical_or(mask_1, mask_2)
    elif fct == "XOR":
        return numpy.logical_xor(mask_1, mask_2)
    else:
        raise ValueError, "function %s is not known"%(fct)
    
def check_bool(array):
    """
    Check if array is boolean and change it if not.
    """

    try:
        mtype = array.dtype
    except:
        __array = numpy.array(array, dtype=bool)
    else:
        if mtype is not bool:
            __array = numpy.array(array, dtype=bool)
        else:
            __array = array
    return array

def generate_array(data, dtype, shape):
    """
        Generates the array corresponding to the data. Data is a list of
        dictionnary. These dictionnaries should, at least, contain the following
        keys : 'pos_x' and 'pos_y' to inform on the position on the array and
        dtype for the value that should come in the array (with some exeptions,
        see parameters below).
        
        Parameters :
            data : list
                List of dictionnaries. All dictionnaries should contain the
                following keys :
                    'pos_x' and 'pos_y' : position in the array
                    ttype : the threshold type
                    dtype : the data to insert in the array
            dtype : str
                The data type to insert in the array. This should be one key
                present in the dictionnaries or :
                'count' : The number of data that are in each array position.
            shape : tuple
                The shape of the returned array.
        
        Returns :
            new_array : 2D numpy.array
                The array generated with the parameters.
        >>> data = [{'pos_x' : 0, 'pos_y' : 0, 'item' : 2},\
                    {'pos_x' : 0, 'pos_y' : 1, 'item' : 2},\
                    {'pos_x' : 0, 'pos_y' : 1, 'item' : 1},\
                    {'pos_x' : 0, 'pos_y' : 2, 'item' : 2},\
                    {'pos_x' : 0, 'pos_y' : 3, 'item' : 2},\
                    {'pos_x' : 1, 'pos_y' : 3, 'item' : 3},\
                    {'pos_x' : 2, 'pos_y' : 0, 'item' : 1},\
                    {'pos_x' : 2, 'pos_y' : 3, 'item' : 2},\
                    ]
        >>> new_data = remove_from_threshold(data, 'item', 1)
        >>> narray = generate_array(new_data, 'count', (4,4))
        >>> print narray
        [[ 1.  2.  1.  1.]
         [ 0.  0.  0.  1.]
         [ 1.  0.  0.  1.]
         [ 0.  0.  0.  0.]]
    """
    new_array = numpy.zeros(shape)
    for item in data:
        res = False
        if dtype == 'count':
            new_array[item['pos_x'], item['pos_y']] += 1
        else:
            # Get the correct value to fill the array
            if dtype in ['fit_length']:
                try:
                    value = item[dtype][0]
                except TypeError:
                    value = None
            else:
                value = item[dtype]
            # If there is no value, nothing to do.
            if value is None:
                res = False
            else:
                res = value > new_array[item['pos_x'], item['pos_y']]
            if res:
                new_array[item['pos_x'], item['pos_y']] = value
    return new_array

def remove_from_threshold(data, atype, avalue):
    """
        Remove from the data the type that are below a value.
        
        >>> import numpy
        >>> data = [{'min' : -0.5, 'max' : 0,'loading_rate' : 0.1, 'dist' : 10 },\
                    {'min' : -0.6, 'max' : -0.1,'loading_rate' : 0.2, 'dist' : 11 },\
                    {'min' : -0.5, 'max' : 0,'loading_rate' : 0.15, 'dist' : 12 },\
                    {'min' : -0.4, 'max' : 0.1,'loading_rate' : 0.12, 'dist' : 9 },\
                    {'min' : -0.3, 'max' : 0.2,'loading_rate' : 0.11, 'dist' : 11 },\
                    {'min' : -0.5, 'max' : 0,'loading_rate' : 0.09, 'dist' : 10 },\
                    {'min' : -0.6, 'max' : -0.1,'loading_rate' : 0.11, 'dist' : 12 },\
                    ]
        >>> new_list = remove_from_threshold(data, 'dist', 10)
        >>> #new_list = remove_from_threshold(data, 'distance', 10)
        >>> new_list = remove_from_threshold([], 'dist', 10)
    """
    if type(avalue) in [float, int]:
        # before, avalue was a single string.
        avalue = [avalue, None]
    if data is None or len(data) == 0:
        return data
    if atype not in data[0].keys():
        raise AttributeError, atype + ' not a valid key for your data.'
    if atype in ['fit_length', 'fit_plength']:
        # In those cases, the data is in a list where the first is the value and
        # the second is the standard deviation like : [value, stdev]
        try:
            if avalue[0] is not None:
                data = [item for item in data if item[atype][0] >= avalue[0]]
        except TypeError:
            pass
        try:
            if avalue[1] is not None:
                data = [item for item in data if item[atype][0] <= avalue[1]]
        except TypeError:
            pass
    else:
        if avalue[0] is not None:
            data = [item for item in data if item[atype] >= avalue[0]]
        if avalue[1] is not None:
            data = [item for item in data if item[atype] <= avalue[1]]
    return data

def flatten(data):
    """
        Flatten an image.
        
        >>> import numpy
        >>> data = numpy.array([range(10) for i in range(10)])
        >>> flatten_data = flatten(data)
        
    """
    vect_x = numpy.arange(data.shape[1])
    # Determine the slope of the drift.
    power = [numpy.polyfit(vect_x, line, 1)[0] for line in data]
    m_pow = numpy.mean(power)
    # Correct the image with the average slope.
    new_data = [line - m_pow * vect_x for line in data]
    new_data = numpy.array(new_data)
    return new_data

def get_message(curr_id, app_name='OpenFovea', app_version='0.1a1'):
    """
        Get message from the server.
    """
    message = None
    headers = {'User-Agent' : '%s/%s'%(app_name, app_version)}
    if curr_id == datetime.datetime(1900, 1, 1):
        # This is the first run.
        # Display then the welcome dialog
        _link = 'http://www.freesbi.ch/docs/openfovea/welcome'
    else:
        _link = 'http://www.freesbi.ch/docs/openfovea/message'
    _req = urllib2.Request(_link, headers=headers)
    try:
        _url = urllib2.urlopen(_req, timeout=0.1)
    except urllib2.URLError:
        _url = None
    if _url is not None:
        _message = _url.readlines()
        if _message[0][:2] == 'id':
            # we had a page, but not the one expected...
            _new_msgid = _message[0].strip('\n').split(' = ')[1]
            if _new_msgid == 'today':
                _new_msgid = datetime.datetime.today()
            else:
                _new_msgid = datetime.datetime.strptime(_new_msgid,
                                                        '%Y, %m, %d')
            if _new_msgid > curr_id:
                message = _message[0:]
                curr_id = _new_msgid
        _url.close()
    return message, curr_id

def ismasked(value):
    '''
        Return TRUE if the value is masked, of FALSE if not.
        Value has to be a scalar.
        
        >>> import numpy
        >>> ismasked(1)
        False
        >>> test = numpy.ma.array([1,2,3,4], mask=[0,0,1,1])
        >>> ismasked(test[0])
        False
        >>> ismasked(test[3])
        True
        
        Only scalar object can be tested with this function :
        >>> ismasked(test)
        Traceback (most recent call last):
        ...
        TypeError: Only scalar object can be tested this way.

    '''
    good_type = True
    try:
        len(value)
        good_type = False
    except:
        pass
    if not good_type:
        raise TypeError, "Only scalar object can be tested this way."
    try:
        result = (type(value) == numpy.ma.core.MaskedConstant)
    except AttributeError:
        result = (type(value) == numpy.ma.core.MaskedArray)
    return result

def export_list(filename, list_array, header, preambule=None, transpose=False):
    """
        Export a list of array in csv.
        
        >>> list = [range(10), range(10), range(20), range(20)]
        >>> header = ['1st', '2nd', '3rd', '4th']
        >>> export_list('test.csv', list, header, transpose=True)
    """
    
    # test that the length of all elements are the same.
    
    length = [len(item) for item in list_array]
    curves = complete_list(list_array)
    if transpose:
        curves = curves.transpose(1,0)
    
    ## generate the filename
    split_filename = os.path.splitext(filename)
    if split_filename[1] in ['.csv', '.txt']:
        filename = split_filename[0]
    else:
        # It is probably not an extention, then restore it
        filename = split_filename[0] + split_filename[1]
    
    filename = filename + '.csv'
    old_masked_print_option = numpy.ma.masked_print_option.display()
    numpy.ma.masked_print_option.set_display('')
    
    # Write the file.
    fid = open(filename, 'wb')
    fcsv = csv.writer(fid)
    if preambule:
        fcsv.writerow(preambule)
    fcsv.writerow(header)
    for item in curves:
        fcsv.writerow(item)
    numpy.ma.masked_print_option.set_display(old_masked_print_option)
    
def complete_list(list_array):
    """
        Complete a list of array in order that all elements has the same
        length.
        
        >>> list = [range(10), range(20)]
        >>> new_list = complete_list(list)
    """
    
    max_length = max([len(item) for item in list_array])
    diff = [max_length - len(item) for item in list_array]
    
    mask = numpy.zeros((len(list_array), max_length))
    array = numpy.zeros((len(list_array), max_length))
    for item, itarray, itmask in zip(list_array, array, mask):
        itarray[:len(item)] = item
        itmask[len(item):] = 1
    new_array = numpy.ma.array(array, mask=mask)
    return new_array
if __name__ == "__main__":
    
    import doctest
    doctest.testmod()
