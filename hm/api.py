#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import math
import string
import datetime
import xarray as xr
import numpy as np
import pandas as pd
import warnings

from .dataarray import HmDomain, HmDataArray, HmSpaceDataArray, HmSpaceTimeDataArray
from .modeltime import ModelTime
from .constants import *
from .utils import *

def set_modeltime(starttime, endtime, timedelta):
    return ModelTime(starttime, endtime, timedelta)    

def set_domain(filename_or_obj, modeltime, varname, is_1d, xy_dimname, **kwargs):
    """Open data array which defines the model domain.
    
    Parameters
    ----------
    filename_or_obj : str 
        String or object which is passed to xarray
    modeltime : HmModelTime
        Object containing the model time
    varname : str, optional
        String specifying the name of the variable which defines
        the model domain
    is_1d : bool
    xy_dimname : str
    
    Returns
    -------
    domain : HmDomain
        The newly created model domain.    
    """
    if is_1d & (xy_dimname is None):
        raise ValueError(
            "If domain dataset is one dimensional then the "
            "name of the space dimension must be specified."
        )
    
    try:
        ds = xr.open_dataset(filename_or_obj, **kwargs)
        if varname is not None:            
            try:
                da = ds[varname]
                # return da

            except KeyError:
                raise KeyError(
                    'File '
                    + filename_or_obj
                    + ' does not contain variable '
                    + varname
                    )

        else:
            if len(ds.data_vars) != 1:
                raise ValueError(
                    'File ' + filename_or_obj + ' contains more than one '
                    'variable. Please specify the variable which '
                    'defines the model domain using the "variable" '
                    'parameter'
                )
            else:
                da = xr.open_dataarray(filename_or_obj, **kwargs)

        if is_temporal(da):
            warnings.warn(
                'Supplied DataArray has a time dimension, but this is '
                'not allowed: taking data from first time point and '
                'discarding the rest.'
            )
            t_dimname = [nm for nm in allowed_t_dim_names if nm in da.dims][0]
            da = da.sel({t_dimname : 0})
            
    except OSError:
        try:
            da = xr.open_rasterio(filename_or_obj, **kwargs)
            return da
        except:
            raise OSError(
                'File ' + filename_or_obj + ' cannot be opened'
            )

    # construct an xarray.Dataset
    dims = get_dimension_names(da, is_1d, xy_dimname)
    coords = get_coordinates(da, dims)
    rename_dict = {value:key for key,value in dims.items()}
    mask = da.astype(bool).rename(rename_dict).rename('mask')
    coords.update({'time' : modeltime.times})
    dims.update({'time' : 'time'})
    ds = xr.merge([xr.Dataset(coords), mask])
    # # TODO: grid cell area should be included in model domain        
    return HmDomain(ds, is_1d=is_1d, xy_dimname=xy_dimname)

def get_files_covering_time_period(filename_or_obj, domain):
    # assume that no more than one file per day is used - is this reasonable?
    datelist = pd.date_range(
        domain.starttime.date(),
        domain.endtime.date(),
        freq='D'
    ).to_pydatetime().tolist()
    filename_list= []
    for date in datelist:
        f = filename_or_obj.format(day=date.day, month=date.month, year=date.year)
        filename_list.append(f)
    # https://stackoverflow.com/a/44628266
    lookup = set()  # a temporary lookup set
    filename_list = [x for x in filename_list if x not in lookup and lookup.add(x) is None]
    return filename_list

def has_format_args(x):
    format_args = [tup[1] for tup in string.Formatter().parse(x) if tup[1] is not None]
    return len(format_args) > 0

def open_xarray_dataset(filename_or_obj, domain, **kwargs):
    if has_format_args(filename_or_obj):
        filename_list = get_files_covering_time_period(filename_or_obj, domain)
    else:
        filename_list = [filename_or_obj]
    kwargs['combine'] = 'by_coords'
    try:
        ds = xr.open_mfdataset(filename_list, **kwargs)        
    except ValueError:
        # try:
        kwargs['decode_times'] = False        
        ds = xr.open_mfdataset(filename_list, **kwargs)
        # does the dataset have a variable which is interpretable as time?            
        timedim = [dim for dim in ds.dims if dim in allowed_t_dim_names][0]
        timevars = [var for var in ds.variables.keys() if var in allowed_t_dim_names and var not in timedim]
        if len(timevars) == 1:
            timevar = timevars[0]
            time = ds[timevar]
            timenum = np.array(
                nc.num2date(time.values, time.units),
                dtype='datetime64'
            )
            ds.update(xr.Dataset({timedim : timenum}))
        else:
            raise ValueError                    
    return ds

def load_hmdataarray(filename_or_obj, **kwargs):
    with open_hmdataarray(filename_or_obj, **kwargs) as da:
        return da.load()

def open_hmdataarray(
        filename_or_obj,
        variable,
        domain,
        is_1d,
        xy_dimname,
        use_xarray=True,
        **kwargs
):
    """Open a dataset from a file or file-like object.

    Parameters
    ----------
    TODO

    Returns
    -------
    TODO

    Notes
    -----
    TODO

    See also
    --------
    TODO
    """
    ds = open_xarray_dataset(filename_or_obj, domain, **kwargs)
    
    da = ds[variable]
    dims = get_dimension_names(da, is_1d, xy_dimname)
    coords = get_coordinates(da, dims)
    if is_1d & (xy_dimname not in da.dims):
        raise ValueError(
            'DataArray object is specified as one-dimensional but does '
            'not contain space dimension ' + xy_dimname
        )
    
    all_dims_in_domain = all(dim in domain.dims for dim in dims.keys())
    if not all_dims_in_domain:
        missing_dims = [value for key,value in dims.items() if key not in domain.dims]
        raise ValueError(
            'DataArray object has dimensions which are not specified in '
            'the model domain: ' + ', '.join(missing_dims)
        )
    
    temporal = is_temporal(da)
    spatial = is_spatial(da, is_1d, xy_dimname)
    if spatial:
        data_is_1d = 'xy' in dims.keys()
        data_is_2d = all(dim in dims.keys() for dim in ('x','y'))
        if data_is_1d:
            if domain.is_2d:
                raise ValueError(
                    'DataArray object is one-dimensional but model domain is '
                    'two-dimensional.'
                )
            
            all_xy = all(xy in coords.xy for xy in domain.coords['xy'])
            if not all_xy:
                raise ValueError(
                    'One-dimensional DataArray does not contain all model '
                    'domain points'
                    )
            
        elif data_is_2d:
            if domain.is_1d:
                raise ValueError(
                    'DataArray object is two-dimensional but model domain is '
                    'one-dimensional.'
                )
                        
            domain_ext = get_extent(domain.coords)
            data_ext = get_extent(da.coords)
            domain_in_data = \
                (data_ext.xmin <= domain_ext.xmin) \
                & (data_ext.xmax >= domain_ext.xmax) \
                & (data_ext.ymin <= domain_ext.ymin) \
                & (data_ext.ymax >= domain_ext.ymax)
            if not domain_in_data:
                raise ValueError(
                    'DataArray does not entirely contain model domain: '
                    'Domain extent is ' + str(tuple(domain_ext.values())) + '(left, right, top, bottom)'
                    'Data extent is ' + str(tuple(data_ext.values())) + '(left, right, top, bottom)'
                )
    if temporal:
        data_starttime = pd.Timestamp(da.coords[dims.time].values[0])
        data_endtime = pd.Timestamp(da.coords[dims.time].values[-1])
        time_domain_in_data = \
            (data_starttime <= domain.starttime) \
            & (data_endtime >= domain.endtime)
        if not time_domain_in_data:
            raise ValueError(
                'DataArray does not entirely contain model time domain: '
                'Domain temporal extent is: ' + str(domain.starttime) + ' -> ' + str(domain.endtime) + ''\
                'Data temporal extent is: ' + str(data_starttime) + ' -> ' + str(data_endtime)
            )

    # rename dimensions to standard names
    rename_dict = {value:key for key,value in dims.items()}
    da = da.rename(rename_dict)
    
    if temporal:
        if spatial:
            hm = HmSpaceTimeDataArray(da, domain, is_1d, xy_dimname)
        else:
            hm = HmTimeDataArray(da, domain, is_1d, xy_dimname)
    else:
        hm = HmSpaceDataArray(da, domain, is_1d, xy_dimname)

    return hm
