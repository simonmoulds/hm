#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import math
import string
import datetime
import xarray as xr
import netCDF4 as nc
import numpy as np
import pandas as pd

# https://github.com/cdgriffith/Box/
from box import Box
from collections import namedtuple, OrderedDict
from .constants import *

def is_temporal(dataarray):
    return np.any([dim in allowed_t_dim_names for dim in dataarray.dims])

def is_spatial(dataarray, is_1d, xy_varname):
    has_x = np.any([nm in allowed_x_dim_names for nm in dataarray.dims])
    has_y = np.any([nm in allowed_y_dim_names for nm in dataarray.dims])
    return (has_x & has_y) | (is_1d & (xy_varname in dataarray.dims))

def get_spatial_extent(coords):
    xres = coords.x[1] - coords.x[0]
    yres = coords.y[1] - coords.y[0]
    xmin = np.min(coords.x) - (xres / 2.)
    xmax = np.max(coords.x) + (xres / 2.)
    ymin = np.min(coords.y) - (yres / 2.)
    ymax = np.max(coords.y) + (yres / 2.)
    return Box(left=xmin, right=xmax, top=ymax, bottom=ymin, frozen_box=True)

def get_dimension_names(dataarray, is_1d=False, xy_dimname=None):
    dimnames = OrderedDict()
    for dim in dataarray.dims:
        if is_1d & (dim == xy_dimname):
            dimnames['xy'] = dim 
        elif dim in allowed_x_dim_names:
            dimnames['x'] = dim
        elif dim in allowed_y_dim_names:
            dimnames['y'] = dim
        elif dim in allowed_z_dim_names:
            dimnames['z'] = dim
        elif dim in allowed_t_dim_names:
            dimnames['time'] = dim
        else:
            dimnames[dim] = dim
    return Box(dimnames, frozen_box=True)

# def get_dimension_axes(dataarray, dimnames):
#     axes = OrderedDict()
#     for dim, dimname in dimnames.items():
#         if dimname is not None:
#             axes[dim] = [position for position,value in enumerate(dataarray.dims) if value == dimname][0]
#     return Box(axes, frozen_box=True)

def get_coordinates(dataarray, dimnames):
    coords = OrderedDict()
    for dim, dimname in dimnames.items():
        if dimname is not None:
            coords[dim] = dataarray[dimname].values
    return Box(coords, frozen_box=True)
    
# NOT USED:
# =========

# def repair_time_string(timestr):    
#     timestr = timestr.lower()
#     timestr_split = timestr.split()
#     units = timestr_split[0]
#     successful = False
#     if timestr_split[1] == 'since':
#         successful = True
#     elif 'since' in timestr_split:
#         time_units = timestr.split('since')[0]
#         allowable_units = ['second','minute','hour','day']
#         occurrences = [unit in time_units for unit in allowable_units]
#         number_of_occurrences = sum(1 for occurrence in occurrences if occurrence is True)
#         if number_of_occurrences == 1:
#             time_unit = [unit for (unit, index) in zip(allowable_units, occurrences) if index][0]
#             timestr = time_unit + ' since' + timestr.split('since')[1]
#             successful = True
#         else:
#             msg = "ambiguous time units"
#     else:
#         msg = "time string does not contain 'since'"
#     if not successful:
#         raise ValueError(msg)    
#     return timestr

# def get_format_args(x):
#     """Function to get format arguments from a string. This is
#     useful to work out whether input data files are provided
#     on a monthly or yearly basis
#     """
#     format_args = [tup[1] for tup in string.Formatter().parse(x) if tup[1] is not None]
#     return format_args

# def check_format_args_ok(format_args, allowable_args, allow_duplicates=True):
#     """Function to check format arguments meet certain
#     conditions.
#     """
#     format_args_ok = all([arg in allowable_args for arg in format_args])
#     if format_args_ok and not allow_duplicates:
#         format_args_ok = not any_duplicates(format_args)
#     return format_args_ok

    
