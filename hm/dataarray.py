#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import math
import string
import datetime
import netCDF4 as nc
import xarray as xr
import numpy as np
import pandas as pd

# https://github.com/cdgriffith/Box/
from box import Box
from collections import namedtuple, OrderedDict

from .constants import *
from .utils import *

# https://gis.stackexchange.com/a/166900


class HmBaseClass(object):
    
    def __init__(
            self,
            dataarray_or_dataset,
            is_1d=False,
            xy_dimname=None,
            model_is_1d=True,
            has_data=True
    ):
        self._data = dataarray_or_dataset
        if is_1d & (xy_dimname not in self._data.dims):
            raise ValueError(
                'DataArray or Dataset is specified as '
                'one-dimensional, but does not contain the '
                'provided space dimension name: ' + xy_dimname
            )
        self._is_1d = is_1d
        self._xy_dimname = xy_dimname
        self._model_is_1d = model_is_1d
        self._has_data = has_data
        self._in_memory = False
        self._update_metadata()

    def _update_metadata(self):
        # Extract metadata from xarray dataset.
        self._dims = get_xr_dimension_names(
            self._data,
            self._is_1d,
            self._xy_dimname
        )
        self._axes = get_xr_dimension_axes(self._data, self._dims)
        self._coords = get_xr_coordinates(self._data, self._dims)
        self._is_1d = ('xy' in self.dims)
        self._is_2d = (
            (not self._is_1d)
            & ('x' in self.dims)
            & ('y' in self.dims)
        )
        self._is_spatial = self._is_1d | self._is_2d
        self._spatial_extent()

    def _spatial_extent(self):
        # Extract spatial extent from data object.

        # Only relevant for two-dimensional datasets; for
        # one-dimensional datasets the extent is None.
        if self.is_2d:
            self._extent = get_spatial_extent(self._coords)
        else:
            self._extent = None

    @property
    def extent(self):
        """Spatial extent of the object."""
        return self._extent

    @property
    def dims(self):
        """Dimensions of the object."""
        return tuple(self._dims.keys())

    @property
    def coords(self):
        """Coordinates of the object."""
        return self._coords

    @property
    def is_1d(self):
        """bool: Whether the object is one-dimensional."""
        return self._is_1d

    @property
    def is_2d(self):
        """bool: Whether or not the object is two-dimensional."""
        return self._is_2d

    @property
    def is_spatial(self):
        """bool: Whether or not the object is spatial."""
        return self.is_1d | self.is_2d

    @property
    def is_temporal(self):
        """bool: Whether or not the object is temporal."""
        return 'time' in self.dims

    @property
    def has_data(self):
        """bool: Whether or not the object has data."""
        return self._has_data

    @property
    def in_memory(self):
        """bool: Whether or not the data is stored in memory."""
        return self._in_memory


class HmDataArray(HmBaseClass):
    pass


class HmDataset(HmBaseClass):
    pass


class HmDomain(HmDataset):
    def __init__(
            self,
            dataarray_or_dataset,
            is_1d=False,
            xy_dimname=None,
            model_is_1d=True,
            has_data=True
    ):
        """Model spatial domain.

        The model spatial domain is the spatial grid over which
        computations are performed.

        To load data, use the ``open_hmdomain`` function.

        Parameters
        ----------
        dataarray_or_dataset: xarray.Dataset or xarray.DataArray
            xarray.Dataset or xarray.DataArray        
        is_1d: bool, optional
            Whether the dataset has a one-dimensional 
            representation of space (i.e. a set of points). This 
            is opposed to a two-dimensional representation which will 
            have coordinates to identify the location of each point.
        xy_dimname: str, optional
            If the dataset is one-dimensional, this 
            parameter specifies the name of the space dimension, which
            is often non-standard (e.g. 'land').
        model_is_1d: bool, optional
            Whether the model is one-dimensional.
        has_data: bool, optional
            Whether or not the xarray object contains data        
        """
        super().__init__(
            dataarray_or_dataset,
            is_1d,
            xy_dimname,
            model_is_1d,
            has_data
        )
        if self._is_2d:
            self._coords['xy'] = np.where(self._data['mask'])

    @property
    def is_latlon(self):
        """bool: Whether or not the domain uses geographical coordinates."""
        return True

    @property
    def y(self):
        """numpy.array: y-coordinates."""
        return self._coords['y']

    @property
    def x(self):
        """numpy.array: x-coordinates."""
        return self._coords['x']

    @property
    def nx(self):
        """int: Size of model grid in x-direction."""
        return len(self._coords['x'])

    @property
    def ny(self):
        """int: Size of model grid in y-direction."""
        return len(self._coords['y'])

    @property
    def nxy(self):
        """int: Total number of model grid points."""
        if self._is_1d:
            return len(self._coords['xy'])
        elif self._is_2d:
            return len(self._coords['xy'][0])

    @property
    def area(self):
        """numpy.array: Area represented by each model grid point."""
        return self._data['area']

    @property
    def mask(self):
        """numpy.array: Model spatial mask.

        This is relevant when the model uses a 2D rectangular grid.
        """
        return self._data['mask']

    # TODO: I don't think these time-related properties should be in HmDomain?
    @property
    def starttime(self):
        """pandas.Timestamp: Model start time."""
        return pd.Timestamp(self._data['time'].values[0])

    @property
    def endtime(self):
        """pandas.Timestamp: Model end time."""
        return pd.Timestamp(self._data['time'].values[-1])

    @property
    def n_timestep(self):
        """int: Number of time points in model simulation."""
        return len(self._coords['time']) - 1

    @property
    def dt(self):
        """pandas.Timedelta: Duration of each timestep."""
        return (self.endtime - self.starttime) / self.n_timestep


class HmSpaceDataArray(HmDataArray):
    def __init__(
            self,
            dataarray,
            domain,
            is_1d=False,
            xy_dimname=None,
            model_is_1d=True,
            has_data=True,
            **kwargs
    ):
        """N-dimensional array with labeled coordinates and dimensions.

        A wrapper for `xarray.DataArray` for datasets which have no time
        dimension, with some additional methods for efficient subsetting 
        based on the model domain. This would typically be used to 
        represent spatially varying input parameters, such as maps of 
        topography, land use/land cover, and soil properties.

        To load data, use the `hm.api.open_hmdataset` function.

        Parameters
        ----------
        dataarray: xarray.DataArray
            xarray.DataArray
        domain: hm.HmDomain
            The model domain.
        is_1d: bool, optional
            Whether the dataset has a one-dimensional 
            representation of space (i.e. a set of points). This 
            is opposed to a two-dimensional representation which will 
            have coordinates to identify the location of each point.
        xy_dimname: str, optional
            If the dataset is one-dimensional, this parameter specifies the name
            of the space dimension, which is often non-standard (e.g. 'land').
        model_is_1d: bool, optional
            Whether the model is one-dimensional.
        has_data: bool, optional
            Whether or not the xarray object contains data        
        **kwargs : Any
            Additional keyword arguments to `xarray.Dataset.sel`
        """
        self._domain = domain
        super().__init__(
            dataarray,
            is_1d,
            xy_dimname,
            model_is_1d,
            has_data
        )
        self.subset(**kwargs)

    def interp(self, **kwargs):
        self._data = self._data.interp(**kwargs)
        self._update_metadata()

    def select(self, **kwargs):
        self._data.sel(**kwargs)
        self._update_metadata()

    def subset(self, skip=[], method='nearest', **kwargs):
        """Subset data with model domain.

        Parameters
        ----------
        skip : list
            List of dimension names to skip when performing subset.
        method : str, optional
            Method to use for inexact matches (see ``xarray.Dataset.sel`` for 
            more details).
        **kwargs : Any
            Additional keyword arguments to `xarray.Dataset.sel`
        """
        self._get_index()
        skip += ['time']
        # skip = skip.append('time')
        # skip.append('time')
        index = {k: v for k, v in self.index.items() if k not in skip}
        if all([dim in self._data.coords for dim in self.index]):
            self._data = self._data.sel(index, method=method, **kwargs)
        else:
            self._data = self._data.sel(index)
        self._update_metadata()

    def _get_index(self):
        # Construct index for subsetting based on domain.
        index_dict = {}
        for dim, dimname in self._dims.items():
            try:
                index_dict[dimname] = self._domain.coords[dim]
            except:
                pass
        self.index = index_dict

    def load(self):        
        """Load data to memory."""
        self._data.load()
        self._in_memory = True

    @property
    def values(self):
        """numpy.array: Array values."""
        # it is a requirement that space dimensions are last - document this (CF Conventions)
        if self._model_is_1d:
            return self._data.values[..., self._domain.mask.values]
        else:
            return self._data.values


# class HmSpaceDataset(object):
#     def __init__(
#             self,
#             dataarray_or_dataset,
#             domain,
#             is_1d=False,
#             xy_dimname=None,
#             model_is_1d=True,
#             has_data=True
#     ):
#         self._domain = domain
#         super().__init__(
#             dataarray_or_dataset,
#             is_1d,
#             xy_dimname,
#             model_is_1d,
#             has_data
#         )
#         self.subset()

class HmSpaceTimeDataArray(HmSpaceDataArray):
    def __init__(
            self,
            dataarray,
            nc_dataset,
            variable,
            domain,
            is_1d=False,
            xy_dimname=None,
            model_is_1d=True,
            has_data=True
    ):
        """N-dimensional array with labeled coordinates and dimensions.

        A wrapper for `xarray.DataArray` for spatiotemporal datasets, with some
        additional methods for efficient subsetting based on the model domain.
        This would typically be used to represent model driving data (e.g. 
        meteorological variables).

        To load data, use the `hm.api.open_hmdataset` function.

        Parameters
        ----------
        dataarray: xarray.DataArray
            xarray.DataArray
        nc_dataset : netCDF4.Dataset
            The netCDF4.Dataset object, providing a lower-level interface to the
            underlying dataset which returns selected data much faster than 
            comparable operations in xarray. The trade-off is that only exact
            indices can be used. To overcome this, hm.HmSpaceTimeDataArray 
            includes a 'ghost' xarray.DataArray object which is used for (exact 
            and inexact) indexing. The output of xarray.DataArray.sel is used to 
            find the exact index in the corresponding netCDF4 object.
        variable : str
            Variable name
        domain: hm.HmDomain
            The model domain.
        is_1d: bool, optional
            Whether the dataset has a one-dimensional 
            representation of space (i.e. a set of points). This 
            is opposed to a two-dimensional representation which will 
            have coordinates to identify the location of each point.
        xy_dimname: str, optional
            If the dataset is one-dimensional, this parameter specifies the name
            of the space dimension, which is often non-standard (e.g. 'land').
        model_is_1d: bool, optional
            Whether the model is one-dimensional.
        has_data: bool, optional
            Whether or not the xarray object contains data        
        """
        super().__init__(
            dataarray,
            domain,
            is_1d,
            xy_dimname,
            model_is_1d,
            has_data
        )
        self._varname = variable
        self._nc_data = nc_dataset  # nc.Dataset(filename_or_obj, 'r')
        self._nc_coords = get_nc_coordinates(self._nc_data, self._dims)
        self._get_nc_index()
        self.select(time=self._domain.starttime, method='nearest')

    def _get_nc_index(self):
        """Map xarray to underlying netCDF4 dataset.

        This is necessary because xarray is rather slow at
        retrieving data from file.
        """
        nc_slice = [slice(None)] * len(self._dims)
        for i, (dim, dimname) in enumerate(self._dims.items()):
            nc_dim = self._nc_coords[dim]
            xr_dim = self.index[dimname]
            index = match(xr_dim, nc_dim)
            nc_slice[i] = index
        self._nc_index = nc_slice

    def select(self, time, **kwargs):
        """Select a temporal subset of the data.

        Parameters
        ----------
        time : pandas.Timestamp, slice
            The time at which data is required. If time is a
            slice, the start and stop points must be
            pandas.Timestamp objects.
        **kwargs : Any
            Additional keyword arguments to
            `xarray.Dataset.sel`
        """
        # select time using xarray selection
        xr_data = self._data.sel({self._dims['time']: time}, **kwargs)
        xr_time = xr_data[self._dims['time']].values

        # match xarray times with those in netCDF4 file (there should always be a match)
        mf_time_index = match(xr_time, self.nc_time)

        # index of times in individual netCDF4 files
        time_index = self._nc_coords['_index'][mf_time_index]
        try:
            time_index.__len__()
        except AttributeError:
            time_index = [time_index]

        # index of files over which the time index is spread
        file_index = np.unique(self._nc_coords['_file'][mf_time_index])
        slc = self._nc_index.copy()

        # 10/2020 - numpy throwing a VisibleDeprecationWarning - 'Creating an ndarray from ragged nested sequences (which is a list-or-tuple of lists-or-tuples-or ndarrays with different lengths or shapes) is deprecated.'
        # converting 'slc' to tuple (i.e. 'tuple(slc)') removes this error
        if len(file_index) == 1:
            slc[self._axes['time']] = time_index
            values = self._nc_data[file_index[0]
                                   ].variables[self._varname][tuple(slc)]
        elif len(file_index) > 1:
            values = []
            for i in file_index:
                slc[self._axes['time']] = time_index[file_index == i]
                values.append(
                    self._nc_data[file_index[i]].variables[self._varname][tuple(slc)])
            values = np.concatenate(values, axis=self._axis['time'])
        # TODO: this throws an UnboundLocalError if values has not yet been created (i.e. empty values)
        self._values = values

    @property
    def values(self):
        """numpy.array: Data values."""
        # TODO: make a separate property - e.g. values_masked
        if self._model_is_1d:
            return self._values[..., self._domain.mask.values]
        else:
            return self._values

    @property
    def nc_time(self):
        """numpy.array: Dataset time coordinates."""
        return self._nc_coords[self._dims['time']]

    @property
    def starttime(self):
        """pandas.Timestamp: Start time of the dataset."""
        return pd.Timestamp(self._data['time'].values[0])

    @property
    def endtime(self):
        """pandas.Timestamp: End time of the dataset."""
        return pd.Timestamp(self._data['time'].values[-1])

    @property
    def n_timestep(self):
        """int: number of time points in the dataset."""
        return len(self._data['time']) - 1

    @property
    def dt(self):
        """pandas.Timedelta: time between each time point."""
        # This assumes a consistent timestep
        return (self.endtime - self.starttime) / self.n_timestep
