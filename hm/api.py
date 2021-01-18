#!/usr/bin/env python
# -*- coding: utf-8 -*-

import string
import xarray as xr
import numpy as np
import pandas as pd
import warnings

warnings.simplefilter(action='ignore', category=FutureWarning)

from .dataarray import HmDomain, HmDataArray, HmSpaceDataArray, HmSpaceTimeDataArray
from .modeltime import ModelTime
from .constants import *
from .utils import *


def set_modeltime(starttime, endtime, timedelta):
    return ModelTime(starttime, endtime, timedelta)

# See the following answer to compute area from lat/long grid:
# https://gis.stackexchange.com/a/232996


class InputError(Exception):
    def __init__(self, expression, message):
        self.expression = expression
        self.message = message

def set_domain(
        filename_or_obj,
        modeltime,
        mask_varname=None,
        area_varname=None,
        is_1d=False,
        xy_dimname=None,
        z_coords=None,
        pseudo_coords=None,
        **kwargs):
    """Define the model domain.

    :param filename_or_obj: String or object which is passed to xarray
    :type filename_or_obj: str
    :param modeltime: Object containing the temporal grid of the model
    :type modeltime: class:`hm.modeltime.HmModelTime`
    :param mask_varname: String specifying the name of the variable 
        which defines the model domain
    :type mask_varname: str, optional
    :param area_varname: String specifying the name of the variable 
        which defines the model domain. If provided and 
        `filename_or_obj` represents a netCDF file, the program 
        searches for the variable in the file, raising an error if it 
        isn't found. If not provided, or if `filename_or_obj` is a 
        GDAL raster file, it is assumed that the values in the mask 
        represent the grid cell area. Note that not all models require 
        the cell area; in these cases it is perfectly fine to use
        an arbitrary value (e.g. 1)
    :type area_varname: str, optional
    :param is_1d: Whether space is represented as a 2-dimensional 
        grid or 1-dimensional set of points
    :type is_1d: bool, optional
    :param xy_dimname: If `is_1d = True`, then `xy_dimname` is the name 
        of the space dimension in `filename_or_obj`
    :type xy_dimname: str, optional
    :param pseudo_coords: Pseudo coordinates are defined as coordinates 
        which are typically used in hydrological models to represent
        sub-grid variables (e.g. land use/land cover types),
        but may be used for other purposes (e.g. ensemble
        simulations). The dictionary keys should match the
        name of the corresponding dimension in any input netCDF
        file, as well as the name used in the `variable_list`
        module for model implementations build on the `hm`
        infrastructure. Dictionary values should be a
        1-dimensional numpy array with the value of each
        coordinate (typically a sequence of integers)
    :type pseudo_coords: dict
    :param **kwargs: Keyword arguments passed to 
        class:`xarray.open_dataset()`
    :type **kwargs: any        

    :return: An object describing the spatiotemporal characteristics 
        of a model domain.
    :rtype: class:`hm.domain.HmDomain`

    """
    if is_1d & (xy_dimname is None):
        raise ValueError(
            'If domain dataset is one dimensional then the '
            'name of the space dimension must be specified.'
        )
    try:
        ds = xr.open_dataset(filename_or_obj, **kwargs)
        if is_temporal(ds):
            warnings.warn(
                'Supplied DataArray has a time dimension, but this is '
                'not allowed: taking data from first time point and '
                'discarding the rest.'
            )
            t_dimname = [nm for nm in allowed_t_dim_names if nm in ds.dims][0]
            ds = ds.isel({t_dimname: 0})

        if mask_varname is not None:
            try:
                mask = ds[mask_varname]
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
                mask = xr.open_dataarray(filename_or_obj, **kwargs)

        mask.name = 'mask'
        if area_varname is not None:
            try:
                area = ds[area_varname]
            except KeyError:
                raise KeyError(
                    'File '
                    + filename_or_obj
                    + ' does not contain variable '
                    + varname
                )
        else:
            # should we include a warning here?
            area = mask.copy()
            area.name = 'area'
            area[:] = 1

        dims = get_xr_dimension_names(mask, is_1d, xy_dimname)
        coords = get_xr_coordinates(mask, dims)

        # If the model is one-dimensional then the domain must
        # contain variables which specify the x and y coordinates
        if is_1d:                                # 1D data has to be specified as netCDF
            varnames = list(ds.variables.keys())
            has_x = np.any([nm in allowed_x_dim_names for nm in varnames])
            has_y = np.any([nm in allowed_y_dim_names for nm in varnames])
            if not (has_x & has_y):
                raise InputError('Input netCDF must have spatial information')            
            xnm = [nm for nm in varnames if nm in allowed_x_dim_names]
            ynm = [nm for nm in varnames if nm in allowed_y_dim_names]
            if (len(xnm) > 1) | (len(ynm) > 1):
                raise InputError('One-dimensional netCDF must have unambiguous variable names for x and y coordinates')
            # Update coords with x and y coordinates
            coords.update(
                {'x' : ds.variables[xnm[0]].values,
                 'y' : ds.variables[ynm[0]].values}
            )            
                 
    except OSError:
        # An OSError would result if filename_or_obj couldn't
        # be opened with xarray.open_dataset(). Next we try to
        # to open the dataset xarray.open_rasterio(), for the
        # case wheree filename_or_obj is a GDAL dataset (must be 2d)
        try:
            mask = xr.open_rasterio(filename_or_obj, **kwargs)
            if 'band' in mask.dims:
                mask = mask.sel(band=1)
            mask.name = 'mask'
            area = mask.copy()
            area.name = 'area'

        except:
            raise OSError(
                'File ' + filename_or_obj + ' cannot be opened'
            )

        dims = get_xr_dimension_names(mask, is_1d, xy_dimname)
        coords = get_xr_coordinates(mask, dims)  # separate 1d/2d methods?
                
    # now we construct an xarray.Dataset to represent the
    # dimensions and pseudo-dimensions
    rename_dict = {value: key for key, value in dims.items()}
    mask = mask.astype(bool).rename(rename_dict)
    area = area.rename(rename_dict)
    coords.update({'time': modeltime.times})
    
    # join z, pseudo coordinates
    if z_coords is None:
        z_coords = {}
    if pseudo_coords is None:
        pseudo_coords = {}

    try:
        all_coords = {**z_coords, **pseudo_coords, **coords}
        # all_coords = {**pseudo_coords, **coords}
    except:
        all_coords = coords
    ds = xr.merge([xr.Dataset(all_coords), mask, area])
    return HmDomain(ds, is_1d=is_1d, xy_dimname='xy')


def _get_files_covering_time_period(filename_or_obj, domain):
    # assume that no more than one file per day is used - is this reasonable?
    datelist = pd.date_range(
        domain.starttime.date(),
        domain.endtime.date(),
        freq='D'
    ).to_pydatetime().tolist()
    filename_list = []
    for date in datelist:
        f = filename_or_obj.format(
            day=date.day, month=date.month, year=date.year)
        filename_list.append(f)
    # https://stackoverflow.com/a/44628266
    lookup = set()  # a temporary lookup set
    filename_list = [
        x for x in filename_list if x not in lookup and lookup.add(x) is None]
    return filename_list


def _get_format_args(x):
    return [tup[1] for tup in string.Formatter().parse(x) if tup[1] is not None]

# def _has_format_args(x):
#     format_args = _get_format_args(x)
#     # format_args = [tup[1]
#     #                for tup in string.Formatter().parse(x) if tup[1] is not None]
#     return len(format_args) > 0

def _has_date_format_args(x):
    format_args = _get_format_args(x)
    return np.any([x in format_args for x in ['year','month','day']])


def _has_mc_format_args(x):
    format_args = _get_format_args(x)
    return np.any([x in format_args for x in ['sample']])


def _get_filename_list(filename_or_obj, domain, sample=1):
    if _has_date_format_args(filename_or_obj):
        filename_list = _get_files_covering_time_period(
            filename_or_obj,
            domain
        )
    elif _has_mc_format_args(filename_or_obj):
        filename_list = [filename_or_obj.format(sample=sample)]
    else:
        filename_list = [filename_or_obj]    

def _open_xarray_dataset(filename_or_obj, domain, **kwargs):
    # if _has_format_args(filename_or_obj):
    #     filename_list = _get_files_covering_time_period(
    #         filename_or_obj, domain)
    # else:
    #     filename_list = [filename_or_obj]
    # if _has_date_format_args(filename_or_obj):
    #     filename_list = _get_files_covering_time_period(
    #         filename_or_obj,
    #         domain
    #     )
    # elif _has_mc_format_args(filename_or_obj):
    #     filename_list = _get_files_for_current_sample(
    #         filename_or_obj,
    #         domain
    #     )
    # else:
    #     filename_list = [filename_or_obj]
    filename_list = _get_filename_list(filename_or_obj, domain)
        
    if len(filename_list) > 1:
        kwargs['combine'] = 'by_coords'
    try:
        ds = xr.open_mfdataset(filename_list, **kwargs)
    except ValueError:
        # A ValueError arises when xarray cannot properly
        # decode the times in the dataset (this happens
        # with WFDEI data, for instance). In this case we
        # have to decode the times manually.
        kwargs['decode_times'] = False
        ds = xr.open_mfdataset(filename_list, **kwargs)
        # does the dataset have a variable which is interpretable as time?
        timedim = [dim for dim in ds.dims if dim in allowed_t_dim_names][0]
        # timevars = [var for var in ds.variables.keys(
        # ) if var in allowed_t_dim_names]# and var not in timedim]
        timevars = [var for var in ds.variables.keys(
        ) if var in allowed_t_dim_names and var not in timedim]
        if len(timevars) == 1:
            timevar = timevars[0]
            time = ds[timevar]
            timenum = np.array(
                nc.num2date(time.values, time.units),
                dtype='datetime64'
            )
            ds.update(xr.Dataset({timedim: timenum}))
        else:
            raise ValueError
    return ds


def _open_netcdf_dataset(filename_or_obj, domain, **kwargs):
    filename_list = _get_filename_list(filename_or_obj, domain)
    # if _has_format_args(filename_or_obj):
    #     filename_list = _get_files_covering_time_period(
    #         filename_or_obj, domain)
    # else:
    #     filename_list = [filename_or_obj]
    
    # MFDataset results in ValueError:
    # ValueError: MFNetCDF4 only works with NETCDF3_* and NETCDF4_CLASSIC formatted files, not NETCDF4.
    # until this is resolved, we have to implement the
    # necessary functionality ourselves
    # ds = nc.MFDataset(filename_list, 'r')
    # N.B. netCDF4 is more robust compared to xarray
    # when it comes to interpreting time
    ds = [nc.Dataset(f, 'r') for f in filename_list]
    return ds


def load_hmdataarray(filename_or_obj, **kwargs):
    with open_hmdataarray(filename_or_obj, **kwargs) as da:
        return da.load()
        
def open_hmdataarray(
        filename_or_obj,
        variable,
        domain,
        is_1d=False,
        xy_dimname=None,
        model_is_1d=True,
        use_xarray=True,
        xarray_kwargs={},
        **kwargs
):
    """Open a dataset from a file or file-like object.

    Parameters
    ----------
    filename_or_obj : str
        String or object which is passed to xarray
    variable : str
        The name of the variable to read
    domain : HmDomain
        The spatio-temporal model domain
    is_1d : bool, optional
        Whether space is represented as a 2-dimensional
        grid or 1-dimensional set of points
    xy_dimname : str, optional
        If `is_1d = True`, then `xy_dimname` is the name of the
        space dimension in `filename_or_obj`
    use_xarray : bool, optional
        Currently not used

    Returns
    -------
    HmDataArray

    Notes
    -----
    TODO

    See also
    --------
    TODO
    """

    # TODO: test whether this is OK:
    ds = _open_xarray_dataset(filename_or_obj, domain, **xarray_kwargs)
    da = ds[variable]
    has_data = True
    # try:
    #     ds = _open_xarray_dataset(filename_or_obj, domain, **xarray_kwargs)
    #     da = ds[variable]
    #     has_data = True
    # except OSError:
    #     da = xr.DataArray(
    #         data=np.full(domain._data['area'].shape, np.nan),
    #         dims=domain._data['area'].dims,
    #         coords=domain._data['area'].coords
    #     )
    #     has_data = False

    dims = get_xr_dimension_names(da, is_1d, xy_dimname)
    coords = get_xr_coordinates(da, dims)
    if is_1d & (xy_dimname not in da.dims):
        raise ValueError(
            'DataArray object is specified as one-dimensional but does '
            'not contain space dimension ' + xy_dimname
        )

    # TODO: evaluate whether we need this check?
    # all_dims_in_domain = all(dim in domain.dims for dim in dims.keys())
    # if not all_dims_in_domain:
    #     missing_dims = [value for key,value in dims.items() if key not in domain.dims]
    #     raise ValueError(
    #         'DataArray object has dimensions which are not specified in '
    #         'the model domain: ' + ', '.join(missing_dims)
    #     )

    temporal = is_temporal(da)
    spatial = is_spatial(da, is_1d, xy_dimname)
    if spatial:
        data_is_1d = 'xy' in dims.keys()
        data_is_2d = all(dim in dims.keys() for dim in ('x', 'y'))
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

            domain_ext = get_spatial_extent(domain.coords)
            data_ext = get_spatial_extent(coords)
            domain_in_data = (data_ext.left <= domain_ext.left) \
                & (data_ext.right >= domain_ext.right) \
                & (data_ext.bottom <= domain_ext.bottom) \
                & (data_ext.top >= domain_ext.top)
            if not domain_in_data:
                raise ValueError(
                    'DataArray does not entirely contain model domain: '
                    'Domain extent is ' +
                    str(tuple(domain_ext.values())) +
                    '(left, right, top, bottom)'
                    'Data extent is ' +
                    str(tuple(data_ext.values())) +
                    '(left, right, top, bottom)'
                )

    else:
        raise InputError('Input netCDF must have spatial information')
    
    if temporal:
        # Check the data covers the time domain
        data_starttime = pd.Timestamp(da.coords[dims.time].values[0])
        data_endtime = pd.Timestamp(da.coords[dims.time].values[-1])
        time_domain_in_data = (data_starttime <= domain.starttime + domain.dt/2) \
            & (data_endtime >= domain.endtime - domain.dt/2)
        if not time_domain_in_data:
            warnings.warn(
                'DataArray from file: ' +
                filename_or_obj + ' ' +
                'does not entirely contain model time domain: '
                'Domain temporal extent is: ' +
                str(domain.starttime) + ' -> ' + str(domain.endtime) + ''
                'Data temporal extent is: ' +
                str(data_starttime) + ' -> ' + str(data_endtime)
            )
            # raise ValueError(
            #     'DataArray does not entirely contain model time domain: '
            #     'Domain temporal extent is: ' +
            #     str(domain.starttime) + ' -> ' + str(domain.endtime) + ''
            #     'Data temporal extent is: ' +
            #     str(data_starttime) + ' -> ' + str(data_endtime)
            # )

    if temporal:
        # If temporal then we also open the netCDF file
        # using the netCDF4 library, which retrieves data
        # from the netCDF file much faster than xarray.
        # The time penalty of using xarray is negligible
        # when a file is only accessed once during the
        # simulation (i.e. a spatial data), but can be
        # considerable when accessed multiple times (i.e.
        # (spatio-)temporal data).
        nc_dataset = _open_netcdf_dataset(filename_or_obj, domain)
        # if spatial:
        hm = HmSpaceTimeDataArray(
            da,
            nc_dataset,
            variable,
            domain,
            is_1d,
            xy_dimname,
            model_is_1d,
            has_data,
            **kwargs
        )
        # else:            
        #     hm = HmTimeDataArray(
        #         da,
        #         nc_dataset,
        #         domain,
        #         has_data,
        #         **kwargs
        #     )
    else:
        hm = HmSpaceDataArray(
            da,
            domain,
            is_1d,
            xy_dimname,
            model_is_1d,
            has_data,
            **kwargs)
    return hm
