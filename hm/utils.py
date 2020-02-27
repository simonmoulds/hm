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
    # TODO: document
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

# def get_dimxension_names(dims, is_1d, xy_dimname):
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
    return Box(coords, frozen_box=False)
    # return Box(coords, frozen_box=True)

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
                time_coord = []
                index_coord = []
                file_coord = []
                for i, _ in enumerate(dataset):
                    # TODO: check assumption that files are ordered by times
                    times = decode_nc_times(dataset[i].variables[dimname])
                    time_coord.append(times)
                    file_coord.append(np.array([i] * len(times), dtype=np.int32))
                    index_coord.append(np.arange(len(times)))
                coords[dim] = np.concatenate(time_coord)
                coords['_file'] = np.concatenate(file_coord)
                coords['_index'] = np.concatenate(index_coord)
            else:
                # TODO: check assumption that files have the same spatial/pseudo coordinates
                coords[dim] = dataset[0].variables[dimname][:].data
                
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

def get_format_args(x):
    """Function to get format arguments from a string. This is
    useful to work out whether input data files are provided
    on a monthly or yearly basis
    """
    format_args = [tup[1] for tup in string.Formatter().parse(x) if tup[1] is not None]
    return format_args

def check_format_args_ok(format_args, allowable_args, allow_duplicates=True):
    """Function to check format arguments meet certain
    conditions.
    """
    format_args_ok = all([arg in allowable_args for arg in format_args])
    if format_args_ok and not allow_duplicates:
        format_args_ok = not any_duplicates(format_args)
    return format_args_ok
    
# contents of old file_handling.py:

# #!/usr/bin/env python
# # -*- coding: utf-8 -*-

# import os
# import shutil
# import subprocess
# import datetime
# import math
# import operator
# import calendar
# import string
# import netCDF4 as nc
# import numpy as np
# import sqlite3
# import rasterio
# import xarray as xr
# import pandas as pd

# from numba import jit
# from functools import reduce

# from .Messages import ModelError, ModelFileError, ModelWarning

# import logging
# logger = logging.getLogger(__name__)

# # file cache to minimize/reduce opening/closing files.  
# file_cache = dict()
# attr_cache = dict()

# # Global variables:
# missing_value = 1e20
# # smallNumber = 1E-39

# # Tuple of netcdf file suffixes (extensions) that can be used:
# netcdf_suffixes = ('.nc4','.nc')

# allowed_y_dim_names = ['lat','latitude','y']
# allowed_x_dim_names = ['lon','longitude','x']
# allowed_t_dim_names = ['time']

# class InputNetCDF(object):
#     def __init__(self, model, nc_filename, nc_varname):
#         self.model = model
#         self.nc_filename = nc_filename
#         self.nc_varname = nc_varname
#         self.open_dataset()
#         self.rename_dimensions()
#         self.subset_dataset()
        
#     def rename_dimensions(self):
#         """Standardize dimension variables."""
#         dim_names = [dim_name for dim in self.dataset.keys()]
#         y_name = [dim for dim in allowed_y_dim_names if dim in dim_names]
#         if len(y_name) == 1:
#             self.dataset.rename({y_name : 'lat'})            
#         x_name = [dim for dim in allowed_x_dim_names if dim in dim_names]
#         if len(x_name) == 1:
#             self.dataset.rename({x_name : 'lon'})
#         t_name = [dim for dim in allowed_x_dim_names if dim in dim_names]
#         if len(t_name) == 1:
#             self.dataset.rename({t_name : 'time'})
        
#     def open_dataset(self):
#         self.dataset = xr.open_dataset(self.nc_filename)[nc_varname]
        
#     def subset_dataset(self):
#         self.dataset_subset = self.dataset.sel(
#             lon=self.model.clone_x_coords,
#             lat=self.model.clone_y_coords,
#             method='nearest'
#         )
        
#     @property
#     def values(self):
#         self.dataset_subset.values
    
# class InputTimeVaryingNetCDF(InputNetCDF):    
#     def open_dataset(self):
#         nc_filename_list = self.get_netcdf_files_covering_time_period()
#         self.dataset = xr.open_mfdataset(nc_filename_list)
        
#     def get_netcdf_files_covering_time_period(self):
#         delta = (self.model._modelTime.endDate - self.model._modelTime.startDate).days + 1
#         datelist = pd.date_range(self.model._modelTime.startDate, periods=delta).to_pydatetime().tolist()
#         filenamelist= []
#         for date in datelist:
#             f = ncFile.format(day=date.day, month=date.month, year=date.year)
#             filenamelist.append(f)
#         # https://stackoverflow.com/a/44628266
#         lookup = set()  # a temporary lookup set
#         filenamelist = [x for x in filenamelist if x not in lookup and lookup.add(x) is None]
#         return filenamelist

#     @property
#     def time_slice(self, start_date, end_date):
#         return self.dataset_subset.sel(time=slice(start_date, end_date))

#     @property
#     def values(self, date):
#         return self.dataset_subset.sel(time=date)

# def clear_cache():
#     """Function to clear file cache."""
#     for item in list(file_cache.keys()):
#         file_cache[item].close()
#         file_cache.pop(item)
#     for item in list(attr_cache.keys()):
#         attr_cache.pop(item)

# def db_connect(db_path):
#     con = sqlite3.connect(db_path)
#     return con

# def get_parameter_id(con, parameter_name):
#     with con:
#         query = "SELECT id FROM parameter_values WHERE name = '" + parameter_name + "'"
#         cur = con.cursor()
#         cur.execute(query)
#         val = cur.fetchone()[0]
#     return val

# def read_crop_parameter_from_sqlite(con, crop_id, parameter_name):
#     parameter_id = get_parameter_id(con, parameter_name)
#     with con:
#         query = "SELECT value,min,max FROM parameter_values WHERE crop_id = " + str(crop_id) + " and parameter_id = " + str(parameter_id)
#         cur = con.cursor()
#         cur.execute(query)
#         val = cur.fetchone()
#     return val

# def read_parameter_from_sqlite(con, parameter_name):
#     parameter_id = get_parameter_id(con, parameter_name)
#     with con:
#         query = "SELECT value,min,max FROM parameter_values WHERE parameter_id = " + str(parameter_id)
#         cur = con.cursor()
#         cur.execute(query)
#         val = cur.fetchone()
#     return val

# def get_format_args(x):
#     """Function to get format arguments from a string. This is
#     useful to work out whether input data files are provided
#     on a monthly or yearly basis
#     """
#     format_args = [tup[1] for tup in string.Formatter().parse(x) if tup[1] is not None]
#     return format_args

# def any_duplicates(lst):
#     """Function to check if a list contains any duplicate
#     items.
#     """
#     return len(lst) != len(set(lst))

# def check_format_args_ok(format_args, allowable_args, allow_duplicates=True):
#     """Function to check format arguments meet certain
#     conditions.
#     """
#     format_args_ok = all([arg in allowable_args for arg in format_args])
#     if format_args_ok and not allow_duplicates:
#         format_args_ok = not any_duplicates(format_args)
#     return format_args_ok

# def read_netCDF(nc_file_name):
#     if nc_file_name in list(file_cache.keys()):
#         f = file_cache[nc_file_name]
#     else:
#         f = nc.Dataset(nc_file_name)
#         f.set_auto_mask(False)
#         file_cache[nc_file_name] = f
#     return f

# def read_GDAL(gdal_file_name, bbox=None):
#     if gdal_file_name in list(file_cache.keys()):
#         f = file_cache[gdal_file_name]
#     else:
#         if bbox != None:
#             gdal_warp(gdal_file_name, 'tmp.tif', bbox[0], bbox[1], bbox[2], bbox[3])
#             f = xr.open_rasterio('tmp.tif')
#             os.remove('tmp.tif')
#         else:
#             f = xr.open_rasterio(gdal_file_name)
#         file_cache[gdal_file_name] = f
#     # f = xr.open_rasterio(gdal_file_name)
#     return f

# def check_if_nc_variable_has_dimension(nc_file_name, variable_name, dimension_name):
#     """Function to check if a variable in a netCDF file has
#     a specified dimension.
#     """
#     result = False
#     try:
#         f = read_netCDF(nc_file_name)
#         result = (dimension_name in f.variables[variable_name].dimensions)
#     except:
#         pass
#     return result

# def check_if_nc_has_variable(nc_file_name, variable_name):
#     logger.debug('Check whether the variable: '+str(variable_name)+' is defined in the file: '+str(nc_file_name))
#     result = False
#     try:
#         f = read_netCDF(nc_file_name)
#         result = str(variable_name) in list(f.variables.keys())
#     except:
#         pass
#     return result

# def get_dimension_variable(nc_file_name, dimension_name):
#     """Function to retrieve a given dimension variable from
#     a netCDF.
#     """
#     if not check_if_nc_has_variable(nc_file_name, dimension_name):
#         dimvar = None
#     else:
#         f = read_netCDF(nc_file_name)
#         dimension_name = str(dimension_name)
#         dimvar = f.variables[dimension_name][:]    
#     return dimvar
    
# def rename_latlong_dims(f, LatitudeLongitude):
#     if LatitudeLongitude == True:
#         try:
#             f.variables['lat'] = f.variables['latitude']
#             f.variables['lon'] = f.variables['longitude']
#         except:
#             pass
#     return f

# def get_time_dimension_name(f):
#     dimnames = list(f.dimensions.keys())
#     time_dimension_name = None
#     if 'time' in dimnames:
#         time_dimension_name = 'time'
#     elif 'tstep' in dimnames:
#         time_dimension_name = 'tstep'
#     else:
#         pass    
#     return time_dimension_name
        
# def get_time_variable_name(f):
#     time_variable_name = None
#     if 'time' in f.variables:
#         time_variable_name = 'time'
#     elif 'timestp' in f.variables:
#         time_variable_name = 'timestp'
#     else:
#         pass    
#     return time_variable_name

# def get_time_units(nctime):
#     return repair_time_string(nctime.units)

# def get_time_calendar(nctime):
#     try:
#         calendar = nctime.calendar
#     except:
#         calendar = 'standard'
#     return calendar
    
# # TODO: refactor; allow 1D arrays
# def resample_nc_data(f, varName, cloneMapFileName, timeDimName = None, timeIndex = None):
#     """Function to compare the spatial attributes of input 
#     map with those of the clone map.
#     """
#     var_dims = f.variables[varName].dimensions
#     slc = [slice(None)] * len(var_dims)
#     if (timeIndex is not None) and (timeDimName is not None):
#         time_axis = [i for i in range(len(var_dims)) if var_dims[i] == timeDimName][0]
#         slc[time_axis] = timeIndex
        
#     cropData = f.variables[varName][slc]
#     input_latitudes = f.variables['lat'][:]
#     input_longitudes = f.variables['lon'][:]  # TODO sort this out - add lat/long to attributes?
#     # TEMPORARY FIX:
#     if len(input_latitudes) > 1:
#         input_attr = get_netcdf_attributes(f)
#         netcdf_y_orientation_follow_cf_convention = (input_latitudes[0] - input_latitudes[1]) > 0
#         if not netcdf_y_orientation_follow_cf_convention:
#             cropData = np.flip(cropData, axis=-2)
#             input_latitudes = input_latitudes[::-1]    
#         if cloneMapFileName != None:
#             clone_map = read_GDAL(cloneMapFileName)
#             clone_attr = get_gdal_attributes(clone_map, cloneMapFileName)
#             sameClone = compare_maps(input_attr, clone_attr)
#             factor = 1
#             if sameClone == False:
#                 logger.debug('Crop to the clone map with upper left corner (x,y): '+ str(clone_attr['xUL']) + ' , ' + str(clone_attr['yUL']))
#                 cropData = crop_data(cropData, input_latitudes, input_longitudes, input_attr, clone_attr)
#                 factor = int(round(float(input_attr['cellsize']) / float(clone_attr['cellsize'])))
#                 cropData = regrid_data_to_finer_grid(factor, cropData, missing_value)
#     return cropData

# def resample_nc_data_1d(f, varName, cloneMapFileName, timeDimName = None, timeIndex = None):
#     pass

# def resample_nc_data_2d(f, varName, cloneMapFileName, timeDimName = None, timeIndex = None):
#     """Function to compare the spatial attributes of input 
#     map with those of the clone map.
#     """
#     var_dims = f.variables[varName].dimensions
#     slc = [slice(None)] * len(var_dims)
#     if (timeIndex is not None) and (timeDimName is not None):
#         time_axis = [i for i in range(len(var_dims)) if var_dims[i] == timeDimName][0]
#         slc[time_axis] = timeIndex
        
#     cropData = f.variables[varName][slc]
#     input_attr = get_netcdf_attributes(f)
#     input_latitudes = f.variables['lat'][:]
#     input_longitudes = f.variables['lon'][:]  # TODO sort this out - add lat/long to attributes?
#     netcdf_y_orientation_follow_cf_convention = (input_latitudes[0] - input_latitudes[1]) > 0
#     if not netcdf_y_orientation_follow_cf_convention:
#         cropData = np.flip(cropData, axis=-2)
#         input_latitudes = input_latitudes[::-1]    
#     if cloneMapFileName != None:
#         clone_map = read_GDAL(cloneMapFileName)
#         clone_attr = get_gdal_attributes(clone_map, cloneMapFileName)
#         sameClone = compare_maps(input_attr, clone_attr)
#         factor = 1
#         if sameClone == False:
#             logger.debug('Crop to the clone map with upper left corner (x,y): '+ str(clone_attr['xUL']) + ' , ' + str(clone_attr['yUL']))
#             cropData = crop_data(cropData, input_latitudes, input_longitudes, input_attr, clone_attr)
#             factor = int(round(float(input_attr['cellsize']) / float(clone_attr['cellsize'])))
#             cropData = regrid_data_to_finer_grid(factor, cropData, missing_value)
#     return cropData

# def crop_data(x, x_lat, x_lon, x_attr, y_attr):    
#     # longitudes are ascending (i.e. W -> E), hence we *add* half
#     # the input cellsize to the clone map western boundary
#     min_x = min(abs(x_lon[:] - (y_attr['xUL'] + 0.5 * x_attr['cellsize'])))
#     xIdxSta = int(np.where(abs(x_lon[:] - (y_attr['xUL'] + 0.5 * x_attr['cellsize'])) == min_x)[0])
#     xIdxEnd = int(math.ceil(xIdxSta + y_attr['cols'] / (x_attr['cellsize'] / y_attr['cellsize'])))        
#     xStep = int((xIdxEnd - xIdxSta) / abs(xIdxEnd - xIdxSta))
#     # latitudes are descending (i.e. N -> S), hence we *subtract*
#     # half the input cellsize from the clone map northern boundary
#     min_y = min(abs(x_lat[:] - (y_attr['yUL'] - 0.5 * x_attr['cellsize'])))
#     yIdxSta = int(np.where(abs(x_lat[:] - (y_attr['yUL'] - 0.5 * x_attr['cellsize'])) == min_y)[0])
#     yIdxEnd = int(math.ceil(yIdxSta + y_attr['rows'] / (x_attr['cellsize'] / y_attr['cellsize'])))
#     yStep = int((yIdxEnd - yIdxSta) / abs(yIdxEnd - yIdxSta))
#     xSlice = slice(xIdxSta, xIdxEnd, xStep)
#     ySlice = slice(yIdxSta, yIdxEnd, yStep)
#     # cropData = cropData[...,ySlice, xSlice]
#     cropData = x[...,ySlice, xSlice]
#     return cropData
    
# def get_message_for_using_nearest_year(ncFile, varName, oldDate, newDate):
#     msg  = "\n"
#     msg += "WARNING related to the netcdf file: "+str(ncFile)+" ; variable: "+str(varName)+" !"+"\n"
#     msg += "The date "+str(oldDate)+" is NOT available. "
#     msg += "The date "+str(newDate.year)+"-"+str(newDate.month)+"-"+str(newDate.day)+" is used."
#     msg += "\n"
#     return msg

# def get_nearest_date_to_year(year, date):
#     """Function to get the year closest to the specified 
#     year if that is not available.
#     """
#     if date.day == 29 and date.month == 2 and calendar.isleap(date.year) and not calendar.isleap(year):
#         date = datetime.datetime(year, date.month, 28)
#     else:
#         date = datetime.datetime(year, date.month, date.day)
#     return date

# def format_date(date, nctime, useDoy):
#     if isinstance(date, str):
#         date = datetime.datetime.strptime(str(date),'%Y-%m-%d')
#     date = datetime.datetime(date.year,date.month,date.day)  # currently only support daily
#     if useDoy == "yearly":
#         date  = datetime.datetime(date.year,int(1),int(1))
#     if useDoy == "monthly":
#         date = datetime.datetime(date.year,date.month,int(1))
#     if useDoy == "yearly" or useDoy == "monthly" or useDoy == "daily_seasonal":
#         first_year_in_nc_file = findFirstYearInNCTime(nctime)
#         last_year_in_nc_file  =  findLastYearInNCTime(nctime)
#         if date.year < first_year_in_nc_file:
#             date = get_nearest_date_to_year(first_year_in_nc_file, date)
#             msg = get_message_for_using_nearest_year(ncFile, varName, dateInput, date)
#             logger.warning(msg)
#         if date.year > last_year_in_nc_file:
#             date = get_nearest_date_to_year(last_year_in_nc_file, date)
#             msg = get_message_for_using_nearest_year(ncFile, varName, dateInput, date)
#             logger.warning(msg)            
#     return date

# def get_date_index_exact(ncFile, date, t_varname, t_unit, t_calendar):
#     datenum = nc.date2num(date, t_unit, t_calendar)
#     timevar = ncFile[t_varname][:]
#     idx = int(np.nonzero(timevar == datenum)[0][0])
#     return idx

# def get_date_index_before(ncFile, date, t_varname, t_unit, t_calendar):
#     datenum = nc.date2num(date, t_unit, t_calendar)
#     timevar = ncFile[t_varname][:]
#     times_before_current_time = np.nonzero(timevar < datenum)
#     idx = int(np.max(times_before_current_time))
#     return idx

# def get_date_index_after(ncFile, date, t_varname, t_unit, t_calendar):
#     datenum = nc.date2num(date, t_unit, t_calendar)
#     timevar = ncFile[t_varname][:]
#     times_after_current_time = np.nonzero(timevar > datenum)
#     idx = int(np.min(times_after_current_time))
#     return idx
        
# def get_date_availability_message(date, available=True):
#     if available:
#         msg = ("The date %s-%s-%s is available. The 'exact' option is \n"
#                "used while selecting netCDF time."
#                % (str(date.year), str(date.month), str(date.day)))
#     else:
#         msg = ("The date %s-%s-%s is NOT available. The 'exact' option \n"
#                "CANNOT be used while selecting netCDF time."
#                % (str(date.year), str(date.month), str(date.day)))
#     return msg

# def get_message_for_using_date_before_or_after(ncFile, varName, date, using_before=True):
#     msg  = "\n"
#     msg += "WARNING related to the netcdf file: "+str(ncFile)+" ; variable: "+str(varName)+" !"+"\n"
#     if using_before:
#         option = 'before'
#     else:
#         option = 'after'
#     msg += "The date "+str(date.year)+"-"+str(date.month)+"-"+str(date.day)+" is NOT available. The '"+option+"' option is used while selecting netcdf time."
#     msg += "\n"
#     return msg
    
# def get_time_index(ncFile, varName, date, t_varname, t_unit, t_calendar):    
#     try:
#         idx = get_date_index_exact(ncFile, date, t_varname, t_unit, t_calendar)
#         msg = get_date_availability_message(date, available=True)
#         logger.debug(msg)
#     except:
#         msg = get_date_availability_message(date, available=False)
#         logger.debug(msg)
#         try:
#             idx = get_date_index_after(ncFile, date, t_varname, t_unit, t_calendar)
#             msg = get_message_for_using_date_before_or_after(ncFile, varName, date, using_before=False)
#         except:
#             idx = get_date_index_before(ncFile, date, t_varname, t_unit, t_calendar)
#             msg = get_message_for_using_date_before_or_after(ncFile, varName, date, using_before=True)
#         # try:
#         #     idx = get_date_index_before(ncFile, date, t_varname, t_unit, t_calendar)
#         #     msg = get_message_for_using_date_before_or_after(ncFile, varName, date, using_before=True)
#         # except:
#         #     idx = get_date_index_after(ncFile, date, t_varname, t_unit, t_calendar)
#         #     msg = get_message_for_using_date_before_or_after(ncFile, varName, date, using_before=False)
#         logger.debug(msg)
#     return idx

# def get_map_clone_error_message(x):
#     msg  = "\n"
#     msg += "ERROR related to the file: "+str(x)+" !"+"\n"
#     msg += "The map does not have the same attributes (extent, "+"\n"
#     msg += "resolution) as the provided model grid"
#     msg += "\n"
    
# # def read_map_clone(x, y, band=1):
# #     # compare = compare_clone(x, y)
# #     x_attr = get_gdal_attributes(x)
# #     y_attr = get_gdal_attributes(y)
# #     identical = compare_maps(x, y)
# #     if identical == True:
# #         ds = xr.open_rasterio(x)  # TODO - this should be in the class method
# #         band_index = np.max([band - 1, 0])
# #         arr = np.array(ds[band_index,...])
# #     else:
# #         msg = get_map_clone_error_message(x)
# #         raise ValueError,msg
# #     return arr

# def gdal_to_array(gdal_file_name, band=1, cloneMapFileName = None):
#     logger.debug('Reading file: ' + str(gdal_file_name))
#     x = read_GDAL(gdal_file_name)
#     if cloneMapFileName != None:
#         y = read_GDAL(cloneMapFileName)        
#         x_attr = get_gdal_attributes(x, gdal_file_name)
#         y_attr = get_gdal_attributes(y, cloneMapFileName)
#         identical = compare_maps(x_attr, y_attr)
#         if not identical:
#             msg = get_map_clone_error_message(x)
#             raise ValueError(msg)
#     arr = np.array(x)
#     index = int(band) - 1
#     arr = arr[index,...]
#     return arr

# def netcdf_to_array(ncFile,
#                     varName,
#                     dateInput,
#                     isOneDimensional = False,
#                     useDoy = None,                    
#                     cloneMapFileName = None,
#                     cloneMapAttributes = None,
#                     LatitudeLongitude = True,
#                     specificFillValue = None):
    
#     if isOneDimensional:
#         arr = netcdf_to_array_1d(
#             ncFile,
#             varName,
#             dateInput,
#             useDoy,                    
#             cloneMapFileName,
#             cloneMapAttributes,
#             LatitudeLongitude,
#             specificFillValue
#         )
#     else:
#         arr = netcdf_to_array_2d(
#             ncFile,
#             varName,
#             dateInput,
#             useDoy,                    
#             cloneMapFileName,
#             cloneMapAttributes,
#             LatitudeLongitude,
#             specificFillValue
#         )
        
#     return arr
            
# def netcdf_to_array_2d(ncFile,
#                        varName,
#                        dateInput,
#                        useDoy,                    
#                        cloneMapFileName,
#                        cloneMapAttributes,
#                        LatitudeLongitude,
#                        specificFillValue):
    
#     logger.debug('Reading variable: ' + str(varName) + ' from the file: ' + str(ncFile))
#     f = read_netCDF(ncFile)
#     varName = str(varName)
#     f = rename_latlong_dims(f, LatitudeLongitude)
#     t_varname = get_time_variable_name(f)
#     t_dimname = get_time_dimension_name(f)
#     date = format_date(dateInput, f.variables[t_varname], useDoy)
#     t_unit = get_time_units(f.variables[t_varname])
#     t_calendar = get_time_calendar(f.variables[t_varname])
#     timeIndex = get_time_index(f, varName, date, t_varname, t_unit, t_calendar)
#     logger.debug('Using date index ' + str(timeIndex))
#     arr = resample_nc_data(f, varName, cloneMapFileName, t_dimname, timeIndex)
#     f = None
#     return arr

# def netcdf_to_array_1d(ncFile,
#                        varName,
#                        dateInput,
#                        useDoy,                    
#                        cloneMapFileName,
#                        cloneMapAttributes,
#                        LatitudeLongitude,
#                        specificFillValue):
#     # clone map coordinates
#     f = xr.open_dataset(ncFile)
#     clone_lats = f['lat']
#     clone_lons = f['lon']
#     logger.debug('Reading variable: ' + str(varName) + ' from the file: ' + str(ncFile))
#     f = xr.open_dataset(ncFile)            
#     pass

# def netcdf_to_arrayWithoutTime(
#         ncFile, varName,
#         cloneMapFileName  = None,
#         LatitudeLongitude = True,
#         specificFillValue = None,
#         absolutePath = None):        
    
#     logger.debug('Reading variable: ' + str(varName) + ' from the file: ' + str(ncFile))
#     f = read_netCDF(ncFile)
#     varName = str(varName)
#     f = rename_latlong_dims(f, LatitudeLongitude)
#     arr = resample_nc_data(f, varName, cloneMapFileName)
#     f = None
#     return arr

# def netcdf_daily_time_slice_to_array(ncFile,
#                                      varName,
#                                      startDate,
#                                      useDoy = None,
#                                      cloneMapFileName = None,
#                                      LatitudeLongitude = True,
#                                      specificFillValue = None):

#     logger.debug('reading variable: ' + str(varName) + ' from the file: ' + str(ncFile))
#     f = read_netCDF(ncFile)
#     varName = str(varName)
#     f = rename_latlong_dims(f, LatitudeLongitude)
#     t_varname = get_time_variable_name(f)
#     t_dimname = get_time_dimension_name(f)
#     t_unit = get_time_units(f.variables[t_varname])
#     t_calendar = get_time_calendar(f.variables[t_varname])
#     startDate = format_date(startDate, f.variables[t_varname], useDoy)
#     endDate = startDate + datetime.timedelta(days=1)
#     lastDateInNC = nc.num2date(f.variables[t_varname][-1], units=t_unit, calendar=t_calendar)
#     startIndex = nc.date2index(datetime.datetime(startDate.year, startDate.month, startDate.day, 0, 0, 0), f.variables[t_varname])
#     if endDate <= lastDateInNC:
#         endIndex = nc.date2index(datetime.datetime(endDate.year, endDate.month, endDate.day, 0, 0, 0), f.variables[t_varname])
#     else:
#         endIndex = f.variables[t_varname].size        
#     timeIndex = np.arange(startIndex, endIndex)
#     arr = resample_nc_data(f, varName, cloneMapFileName, t_dimname, timeIndex)
#     f = None
#     return arr

# def get_netcdf_files_covering_time_period(ncFile, startDate, endDate):
#     delta = (endDate - startDate).days + 1
#     datelist = pd.date_range(startDate, periods=delta).to_pydatetime().tolist()
#     filenamelist= []
#     for date in datelist:
#         f = ncFile.format(day=date.day, month=date.month, year=date.year)
#         filenamelist.append(f)
#     # https://stackoverflow.com/a/44628266
#     lookup = set()  # a temporary lookup set
#     filenamelist = [x for x in filenamelist if x not in lookup and lookup.add(x) is None]
#     return filenamelist
    
# def netcdf_time_slice_to_array(ncFileUnformatted,
#                                varName,
#                                startDate,
#                                endDate,
#                                useDoy = None,
#                                cloneMapFileName = None,
#                                LatitudeLongitude = True,
#                                specificFillValue = None):

#     # get the list of netcdf files which cover the time period specified by startDate and endDate
#     ncFileList = get_netcdf_files_covering_time_period(
#         ncFileUnformatted, startDate, endDate
#     )

#     # loop though netcdf files
#     arr_list = []
#     for ncFile in ncFileList:        
#         logger.debug('reading variable: ' + str(varName) + ' from the file: ' + str(ncFile))
#         f = read_netCDF(ncFile)
#         varName = str(varName)
#         f = rename_latlong_dims(f, LatitudeLongitude)
#         t_varname = get_time_variable_name(f)
#         t_dimname = get_time_dimension_name(f)
#         t_unit = get_time_units(f.variables[t_varname])
#         t_calendar = get_time_calendar(f.variables[t_varname])
#         startDate = format_date(startDate, f.variables[t_varname], useDoy)
#         endDate = format_date(endDate, f.variables[t_varname], useDoy)
#         lastDateInNC = nc.num2date(f.variables[t_varname][-1], units=t_unit, calendar=t_calendar)
#         startIndex = nc.date2index(
#             datetime.datetime(
#                 startDate.year,
#                 startDate.month,
#                 startDate.day),
#             f.variables[t_varname],
#             select='after'
#         )
#         if endDate <= lastDateInNC:
#             endIndex = nc.date2index(datetime.datetime(endDate.year, endDate.month, endDate.day, 0, 0, 0), f.variables[t_varname], select='after')
#         else:
#             endIndex = f.variables[t_varname].size
#             endIndex -= 1       # account for zero indexing
        
#         timeIndex = np.arange(startIndex, endIndex + 1)
#         arr = resample_nc_data(f, varName, cloneMapFileName, t_dimname, timeIndex)
#         arr_list.append(arr)
#         f = None
        
#     arr = np.concatenate(arr_list, axis=0)
#     return arr
        
# # def compare_clone(inputMapFileName,cloneMapFileName):    
# #     input_attributes = get_map_attributes(inputMapFileName)
# #     clone_attributes = get_map_attributes(cloneMapFileName)
# #     same_clone = True
# #     if input_attributes['cellsize'] != clone_attributes['cellsize']: same_clone = False
# #     if input_attributes['rows'] != clone_attributes['rows']: same_clone = False
# #     if input_attributes['cols'] != clone_attributes['cols']: same_clone = False
# #     if input_attributes['xUL'] != clone_attributes['xUL']: same_clone = False
# #     if input_attributes['yUL'] != clone_attributes['yUL']: same_clone = False    
# #     return same_clone

# def compare_maps(x_attr, y_attr):
#     identical = True
#     if x_attr['cellsize'] != y_attr['cellsize']: identical = False
#     if x_attr['rows'] != y_attr['rows']: identical = False
#     if x_attr['cols'] != y_attr['cols']: identical = False
#     if x_attr['xUL'] != y_attr['xUL']: identical = False
#     if x_attr['yUL'] != y_attr['yUL']: identical = False    
#     return identical

# # def compare_netcdf_with_gdal(nc_file_name, gdal_file_name):
# #     nc_attr = get_netcdf_attributes(nc_file_name)
# #     gdal_attr = get_gdal_attributes(gdal_file_name)
# #     identical = compare_maps(nc_attr, gdal_attr)
# #     return identical

# # def compare_netcdf_with_netcdf(nc1_file_name, nc2_file_name):
# #     nc1_attr = get_netcdf_attributes(nc1_file_name)
# #     nc2_attr = get_netcdf_attributes(nc2_file_name)
# #     identical = compare_maps(nc1_attr, nc2_attr)
# #     return identical

# # def compare_gdal_with_gdal(gdal1_file_name, gdal2_file_name):
# #     gdal1_attr = get_gdal_attributes(gdal1_file_name)
# #     gdal2_attr = get_gdal_attributes(gdal2_file_name)
# #     identical = compare_maps(gdal1_attr, gdal2_attr)
# #     return identical

# def get_gdal_attributes(gdal_file, gdal_file_name):
#     if gdal_file_name in list(attr_cache.keys()):
#         map_attributes = attr_cache[gdal_file_name]
#     else:
#         cellsize = gdal_file.res[0]
#         nrows = gdal_file.y.size
#         ncols = gdal_file.x.size
#         x_upper_left = gdal_file.transform[2]
#         y_upper_left = gdal_file.transform[5]    
#         # if arc_degree == True:
#         #     cellsize = round(cellsize * 360000.)/360000.    
#         map_attributes = {
#             'cellsize' : float(cellsize),
#             'rows'     : float(nrows),
#             'cols'     : float(ncols),
#             'xUL'      : float(x_upper_left),
#             'yUL'      : float(y_upper_left)
#         }
#         attr_cache[gdal_file_name] = map_attributes
        
#     return map_attributes

# def get_gdal_extent(gdal_file, gdal_file_name):
#     map_attributes = get_gdal_attributes(gdal_file, gdal_file_name)
#     xmin = map_attributes['xUL']
#     ymin = map_attributes['yUL'] - map_attributes['rows'] * map_attributes['cellsize']
#     xmax = map_attributes['xUL'] + map_attributes['cols'] * map_attributes['cellsize']
#     ymax = map_attributes['yUL']
#     return (xmin, ymin, xmax, ymax)

# def gdal_warp(input_filename, output_filename, minx, miny, maxx, maxy):
#     """This function is necessary because the python gdal bindings seem
#     unreliable
#     """
#     cmd = ('gdalwarp -overwrite -te '
#            + str(minx) + ' ' + str(miny) + ' ' + str(maxx) + ' ' + str(maxy) + ' '
#            + str(input_filename) + ' '
#            + str(output_filename))
#     subprocess.check_output(cmd, shell=True)
#     # os.system(cmd)

# # def get_netcdf_attributes(nc_file_name):
# #     f = read_netCDF(nc_file_name)
# def get_netcdf_attributes(nc_file):
#     input_latitudes = nc_file.variables['lat'][:]
#     input_longitudes = nc_file.variables['lon'][:]
#     cellsize = float(abs(input_latitudes[0] - input_latitudes[1]))
#     map_attributes = {
#         'cellsize' : cellsize,
#         'rows'     : len(input_latitudes),
#         'cols'     : len(input_longitudes),
#         'xUL'      : np.min(input_longitudes) - 0.5 * cellsize,
#         'yUL'      : np.max(input_latitudes) + 0.5 * cellsize
#     }
#     return map_attributes

# # def get_map_attributes(x):
# #     if isinstance(x, str):
# #         pass
# #     return(x)

# # def get_map_attributes(clone_map, arc_degree=True):
# #     ds = xr.open_rasterio(clone_map)
# #     cellsize = ds.res[0]
# #     nrows = ds.y.size
# #     ncols = ds.x.size
# #     x_upper_left = ds.transform[2]
# #     y_upper_left = ds.transform[5]    
# #     if arc_degree == True:
# #         cellsize = round(cellsize * 360000.)/360000.    
# #     map_attributes = {
# #         'cellsize' : float(cellsize),
# #         'rows'     : float(nrows),
# #         'cols'     : float(ncols),
# #         'xUL'      : float(x_upper_left),
# #         'yUL'      : float(y_upper_left)
# #     }
# #     return map_attributes

# # def getLastDayOfMonth(date):
# #     ''' returns the last day of the month for a given date '''
# #     if date.month == 12:
# #         return date.replace(day=31)
# #     return date.replace(month=date.month + 1, day=1) - datetime.timedelta(days=1)

# def regrid_data_to_finer_grid(factor, coarse, missing_value):
#     # TODO: check with >3 dimensions (only checked with 2D
#     # and 3D arrays so far)    
#     if factor == 1: return coarse
#     dim = np.shape(coarse)
#     ndim = len(dim)
#     nr,nc = dim[-2:]
#     if ndim > 2:
#         othdim = dim[0:(ndim-2)]
#         othdimprod = reduce(operator.mul, othdim)
#         fine = np.zeros(nr * nc * factor * factor * othdimprod).reshape(othdim + (nr * factor, nc * factor)) + missing_value
#     else:
#         fine = np.zeros(nr * nc * factor * factor).reshape(nr * factor, nc * factor) + missing_value

#     dimF = np.shape(fine)
#     nrF,ncF = dimF[-2:]
#     ii = -1
#     for i in range(0, nrF):
#         if i % factor == 0: ii += 1
#         # repeat along space axis, which (I think!!!) should
#         # always be the last dimension of coarse[...,ii,:]
#         # hence, ndim - 2 (-1 to comply with Python zero
#         # indexing, -1 because we select a single row and hence
#         # reduce the number of dimensions by 1)
#         fine[...,i,:] = coarse[...,ii,:].repeat(factor, ndim - 2)

#     return fine

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
            
# def findLastYearInNCTime(ncTimeVariable):    
#     last_datetime = nc.num2date(
#         ncTimeVariable[len(ncTimeVariable) - 1],
#         repair_time_string(ncTimeVariable.units),
#         # ncTimeVariable.units,
#         ncTimeVariable.calendar) 
#     return last_datetime.year

# def findFirstYearInNCTime(ncTimeVariable):
#     first_datetime = nc.num2date(
#         ncTimeVariable[0],
#         repair_time_string(ncTimeVariable.units),
#         ncTimeVariable.calendar)    
#     return first_datetime.year
