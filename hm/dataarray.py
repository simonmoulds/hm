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
    """Base class for the hm data classes."""
    def __init__(
            self,
            dataarray_or_dataset,
            is_1d=False,
            xy_dimname=None
    ):
        """To load data, use the `open_hmdataset` function.

        Parameters:
        -----------
        dataarray_or_dataset : xarray Dataset or DataArray
        is_1d : bool
            Whether the dataset has a one-dimensional representation
            of space (i.e. a set of points). This is opposed to a 
            two-dimensional representation which will have coordinates 
            to identify the location of each point.
        xy_dimname : str 
            If the dataset is one-dimensional, this parameter 
            specifies the name of the space dimension, which is 
            often non-standard (e.g. 'land').
        """
        self._data = dataarray_or_dataset
        if is_1d & (xy_dimname not in self._data.dims):
            raise ValueError(
                'DataArray or Dataset is specified as '
                'one-dimensional, but does not contain the '
                'provided space dimension name: ' + xy_dimname
            )
        self._is_1d = is_1d
        self._xy_dimname = xy_dimname
        self._in_memory = False
        self._update_metadata()
        
    def _update_metadata(self):
        """Extract metadata from xarray dataset."""
        self._dims = get_xr_dimension_names(self._data, self._is_1d, self._xy_dimname)
        self._axes = get_xr_dimension_axes(self._data, self._dims)
        self._coords = get_xr_coordinates(self._data, self._dims)
        self._is_1d = ('xy' in self.dims)
        self._is_2d = (not self._is_1d) & all(dim in self.dims for dim in (self._dims['x'], self._dims['y']))
        self._is_spatial = self._is_1d | self._is_2d
        self._spatial_extent()
        
    def _spatial_extent(self):
        """Extract spatial extent from data object.
        
        Only relevant for two-dimensional datasets; for 
        one-dimensional datasets the extent is None.
        """
        if self.is_2d:
            self._extent = get_spatial_extent(self._coords)
        else:
            self._extent = None            

    @property
    def extent(self):
        return self._extent
    @property
    def dims(self):
        return tuple(self._dims.keys())
    @property
    def coords(self):
        return self._coords
    @property
    def is_1d(self):
        return self._is_1d
    @property
    def is_2d(self):
        return self._is_2d
    @property
    def is_spatial(self):
        return self.is_1d | self.is_2d
    @property
    def is_temporal(self):
        return 'time' in self.dims
    @property
    def in_memory(self):
        return self._in_memory
    
class HmDataArray(HmBaseClass):
    pass

class HmDataset(HmBaseClass):
    pass
        
class HmDomain(HmDataset):
    """Class to represent a model domain.
    
    The model domain is the spatial and temporal grid over 
    which computations are performed.
    """
    @property
    def area(self):
        return self._data['area']
    @property
    def mask(self):
        return self._data['mask']
    @property
    def starttime(self):
        return pd.Timestamp(self._data['time'].values[0])
    @property
    def endtime(self):
        return pd.Timestamp(self._data['time'].values[-1])
    @property
    def nx(self):
        return len(self._coords['x'])
    @property
    def ny(self):
        return len(self._coords['y'])
    @property
    def nxy(self):
        return len(self._coords['xy'])
    @property
    def n_timestep(self):
        return len(self._coords['time']) - 1
    @property
    def dt(self):
        return (self.endtime - self.starttime) / self.n_timestep
    
class HmSpaceDataArray(HmDataArray):
    """Class to represent spatial data."""
    def __init__(
            self,
            dataarray,
            domain,
            is_1d=False,
            xy_dimname=None
    ):
        """To load data, use the `open_hmdataset` function.

        Parameters:
        -----------
        domain : HmDomain object        
        """
        self._domain = domain
        super().__init__(
            dataarray,
            is_1d,
            xy_dimname
        )        
        self.subset()
        
    def subset(self, method='nearest', **kwargs):
        """Subset data with model domain.
        
        Parameters
        ----------
        method : str 
            Method to use for inexact matches (see 
            `xarray.Dataset.sel` for more details).
        **kwargs : Any
            Additional keyword arguments to 
            `xarray.Dataset.sel`
        """
        self._get_index()
        if all([dim in self._data.coords for dim in self.index]):
            self._data = self._data.sel(self.index, method=method, **kwargs)
        else:
            self._data = self._data.sel(self.index)            
        self._update_metadata()
        
    def _get_index(self):
        """Construct index for subsetting based on domain."""
        index_dict = {}
        for dim, dimname in self._dims.items():
            index_dict[dimname] = self._domain.coords[dim]
        self.index = index_dict
        
    def select(self):
        pass
    
    def load(self):
        """Load data to memory."""
        self._data.load()
        self._in_memory = True

class HmSpaceTimeDataArray(HmSpaceDataArray):
    """Class to represent spatio-temporal data."""
    def __init__(
            self,
            dataarray,
            nc_dataset,
            variable,
            domain,
            is_1d,
            xy_dimname
    ):
        """To load data, use the `open_hmdataset` function.

        Parameters:
        -----------
        nc_dataset : netCDF4.Dataset
            The netCDF4.Dataset object, providing a 
            lower-level interface to the dataset which 
            returns selected data much faster than comparable 
            operations in xarray. The trade-off is that only 
            exact indices can be used. To overcome this, this 
            hm.HmSpaceTimeDataArray includes a 'ghost' 
            xarray.DataArray object which is for (exact and
            inexact) indexing. The output of 
            xarray.DataArray.sel is used to find the exact 
            index in the corresponding netCDF4 object. 
        variable : str 
            Variable name
        """
        super().__init__(dataarray, domain, is_1d, xy_dimname)
        self._varname = variable
        self._nc_data = nc_dataset# nc.Dataset(filename_or_obj, 'r')
        self._nc_coords = get_nc_coordinates(self._nc_data, self._dims)
        self._get_nc_index()
        
    def _get_nc_index(self):
        """Map xarray to underlying netCDF4 dataset.
        
        This is necessary because xarray is rather slow at 
        retrieving data from file.
        """
        nc_slice = [slice(None)] * len(self._dims)
        for i, (dim, dimname) in enumerate(self._dims.items()):
            nc_dim = self._nc_coords[dim]
            xr_dim = self.index[dimname]
            nc_slice[i] = match(xr_dim, nc_dim)
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
        xr_data = self._data.sel({self._dims['time'] : time}, **kwargs)
        xr_time = xr_data[self._dims['time']].values
        slc = self._nc_index
        # match xarray times with those in netCDF4 file (there should always be a match)
        slc[self._axes['time']] = match(xr_time, self.nc_time)
        return self._nc_data.variables[self._varname][slc]

    def to_netcdf(self, path, **kwargs):
        # TODO: wrapper to open netcdf file
        try:
            ds = nc.Dataset(path, 'r+')
        except FileNotFoundError:
            ds = nc.Dataset(path, 'w')
        pass
    
    @property
    def nc_time(self):
        return self._nc_coords[self._dims['time']]
    @property
    def starttime(self):
        return pd.Timestamp(self._data['time'].values[0])
    @property
    def endtime(self):
        return pd.Timestamp(self._data['time'].values[-1])
    @property
    def n_timestep(self):
        return len(self._data['time']) - 1
    @property
    def dt(self):
        return (self.endtime - self.starttime) / self.n_timestep
