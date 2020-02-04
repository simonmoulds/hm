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

# file cache for reading netCDF files (this prevents files
# from being reopened at every timestep)
file_cache = dict()

def clear_cache():
    """Function to clear file cache."""
    for item in list(file_cache.keys()):
        file_cache[item].close()
        file_cache.pop(item)

def open_netcdf(filename, **kwargs):
    if filename in list(file_cache.keys()):
        f = file_cache[filename]
    else:
        f = nc.Dataset(filename, **kwargs)
        # f.set_auto_mask(False)
        file_cache[filename] = f
    return f
        
# def add_time_dimension(netcdf, dimname, dimvar, **kwargs):
#     """Add time dimension to a netCDF file.

#     Parameters
#     ----------
#     netcdf : netCDF4.Dataset
#     dimname : str    
#     """
#     shortname = variable_list.netcdf_short_name[dimname]
#     try:
#         datatype = variable_list.netcdf_datatype[dimname]
#     except:
#         datatype = 'f4'
#     dimensions = variable_list.netcdf_dimensions[dimname]
#     netcdf.createDimension(shortname, None)
#     var = netcdf.createVariable(
#         shortname,
#         datatype,
#         dimensions,
#         **kwargs
#     )
#     var.standard_name = variable_list.netcdf_standard_name[dimname]
#     var.long_name = variable_list.netcdf_long_name[dimname]
#     var.units = variable_list.netcdf_unit[dimname]
#     var.calendar = variable_list.netcdf_calendar[dimname]

# def add_nontime_dimension(netcdf, dimname, dimvar, **kwargs):
#     """Add a space or pseudo dimension to a netCDF file.

#     Parameters
#     ----------
#     netcdf : netCDF4.Dataset
#     dimname : str 
#     dimvar : numpy.array 
#     """
#     ndim = len(dimvar)
#     shortname = variable_list.netcdf_short_name[dimname]
#     try:
#         datatype = variable_list.netcdf_datatype[dimname]
#     except:
#         datatype = 'f4'
#     dimensions = variable_list.netcdf_dimensions[dimname]
#     standard_name = variable_list.netcdf_standard_name[dimname]
#     # ensure appropriate number of decimal places
#     # are used for latitude/longitude dimensions
#     if standard_name in ['latitude','longitude']:
#         kwargs['zlib'] = True
#         kwargs['least_significant_digit'] = 16

#     # create the dimension
#     netcdf.createDimension(shortname, ndim)
#     var = netcdf.createVariable(
#         shortname,
#         datatype,
#         dimensions,
#         **kwargs)
#     var.standard_name = variable_list.netcdf_standard_name[dimname]
#     var.long_name = variable_list.netcdf_long_name[dimname]
#     var.units = variable_list.netcdf_unit[dimname]
#     var[:] = np.array(dimvar)
    
def get_variable_dimensions(varname):
    """Get the dimensions of a variable."""
    if isinstance(varname, str):
        var_dims = self.variable_list.netcdf_dimensions[varname]
    elif isinstance(varname, list):
        var_dims = []
        for item in varname:
            var_dims += list(self.variable_list.netcdf_dimensions[item])
        var_dims = tuple(set(var_dims))            
    return var_dims

def reshape_array(arr, mask):
    """Reshape a one-dimensional array.

    This function allocates the values of an array with one 
    space dimension to the corresponding grid cells of a 
    one- or two-dimensional mask.

    Parameters
    ----------
    arr : numpy.array 
    mask : numpy.array, boolean
    """
    arr1 = np.zeros(arr.shape[:-1] + (mask.shape))
    arr1[..., mask] = arr
    return np.ma.array(arr, mask=np.broadcast_to(mask, arr.shape))

def match(x, table):
    table_sorted = np.argsort(table)
    x_pos = np.searchsorted(table[table_sorted], x)
    return table_sorted[x_pos]

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

# def get_dimension_names(dims, is_1d, xy_dimname):
#     dimnames = OrderedDict()
#     for dim in dims:
#         if is_1d & (dim == xy_dimname):
#             dimnames['xy'] = dim 
#         elif dim in allowed_x_dim_names:
#             dimnames['x'] = dim
#         elif dim in allowed_y_dim_names:
#             dimnames['y'] = dim
#         elif dim in allowed_z_dim_names:
#             dimnames['z'] = dim
#         elif dim in allowed_t_dim_names:
#             dimnames['time'] = dim
#         else:
#             dimnames[dim] = dim
#     return Box(dimnames, frozen_box=True)
    
# def get_nc_dimension_names(dataset, is_1d=False, xy_dimname=None):
#     dims = dataset.dimensions
#     return get_dimension_names(dims, is_1d, xy_dimname)

# def get_xr_dimension_names(dataset_or_dataarray, is_1d=False, xy_dimname=None):
#     dims = dataset_or_dataarray.dims
#     return get_dimension_names(dims, is_1d, xy_dimname)
    
def get_xr_dimension_names(dataset_or_dataarray, is_1d=False, xy_dimname=None):
    dimnames = OrderedDict()
    for dim in dataset_or_dataarray.dims:
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

def get_xr_dimension_axes(dataset_or_dataarray, dimnames):
    axes = OrderedDict()
    for dim, dimname in dimnames.items():
        if dimname is not None:
            axes[dim] = [position for position,value in enumerate(dataset_or_dataarray.dims) if value == dimname][0]
    return Box(axes, frozen_box=True)

def get_xr_coordinates(dataset_or_dataarray, dimnames):
    coords = OrderedDict()
    for dim, dimname in dimnames.items():
        if dimname is not None:
            coords[dim] = dataset_or_dataarray[dimname].values
    return Box(coords, frozen_box=True)

def decode_nc_times(timevar):
    try:
        calendar = timevar.calendar
    except AttributeError:
        calendar = 'standard'            
    times = nc.num2date(timevar[:], timevar.units, calendar)
    return np.array(times, dtype='datetime64')

def get_nc_coordinates(dataset, dimnames):
    coords = OrderedDict()
    for dim, dimname in dimnames.items():
        if dimname is not None:
            if dim is 'time':
                coords[dim] = decode_nc_times(dataset.variables[dimname])
            else:
                coords[dim] = dataset.variables[dimname][:].data
                
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

    
