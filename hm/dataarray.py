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

# TODO:
# https://gis.stackexchange.com/a/166900
class HmBaseClass(object):
    def __init__(self, dataarray_or_dataset, is_1d=False, xy_dimname=None):
        self._data = dataarray_or_dataset
        if is_1d & (xy_dimname not in self._data.dims):
            raise ValueError(
                'DataArray or Dataset does not contain the '
                'dimension name: ' + xy_dimname
            )
        self._is_1d = is_1d
        self._xy_dimname = xy_dimname
        self._in_memory = False
        # if self._is_1d:
        #     self._data.rename({xy_dimname : 'xy'})            
        self._update()
        
    def _update(self):
        self._dims = get_xr_dimension_names(self._data, self._is_1d, self._xy_dimname)
        # self._reorder_dims()
        self._axes = get_xr_dimension_axes(self._data, self._dims)
        self._coords = get_xr_coordinates(self._data, self._dims)
        self._is_1d = ('xy' in self.dims)
        self._is_2d = (not self._is_1d) & all(dim in self.dims for dim in (self._dims['x'], self._dims['y']))
        # self._is_2d = (not self._is_1d) & all(dim in self.dims for dim in ('x','y'))
        self._is_spatial = self._is_1d | self._is_2d
        self._spatial_extent()
        
    # def _reorder_dims(self):
    #     pass
    def _spatial_extent(self):
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
    def n_timestep(self):
        return len(self._data['time']) - 1
    @property
    def dt(self):
        return (self.endtime - self.starttime) / self.n_timestep
    
class HmSpaceDataArray(HmDataArray):
    def __init__(self, dataarray, domain, is_1d=False, xy_dimname=None):
        self._domain = domain
        super().__init__(
            dataarray,
            is_1d,
            xy_dimname
        )        
        self.subset()        
    # def _reorder_dims(self):
    #     """Reorder dimensions to match model domain."""
    #     neworder = [dim for dim in self._domain.dims if dim in self.dims]
    #     self._data = self._data.transpose(*neworder)        
    def subset(self):
        """Subset data with model domain."""
        self._get_index()
        # TODO: think about how we can adjust the method arg to sel should vary?
        # TODO: subset in stages; first subset dims with coordinates (so that
        #       a 'method' argument can be supplied, then subset dims which do
        #       not have coordinates.
        if all([dim in self._data.coords for dim in self.index]):
            self._data = self._data.sel(self.index, method='nearest')
        else:
            self._data = self._data.sel(self.index)            
        self._update()
        
    def _get_index(self):
        index_dict = {}
        for dim, dimname in self._dims.items():
            index_dict[dimname] = self._domain.coords[dim]
        self.index = index_dict
        
    def select(self):
        pass
    
    def load(self):
        self._data.load()
        self._in_memory = True

def match(x, table):
    table_sorted = np.argsort(table)
    x_pos = np.searchsorted(table[table_sorted], x)
    return table_sorted[x_pos]

class HmSpaceTimeDataArray(HmSpaceDataArray):
    def __init__(self, dataarray, filename_or_obj, variable, domain, is_1d, xy_dimname):
        super().__init__(dataarray, domain, is_1d, xy_dimname)
        self._varname = variable
        self._nc_data = nc.Dataset(filename_or_obj, 'r')
        self._nc_coords = get_nc_coordinates(self._nc_data, self._dims)
        self._get_nc_index()
        self._nc_time = self._nc_coords[self._dims['time']]
        
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
        time : slice, ...            
        """
        # select time using xarray selection
        xr_data = self._data.sel({self._dims['time'] : time}, **kwargs)
        xr_time = xr_data[self._dims['time']].values
        slc = self._nc_index
        # match xarray times with those in netCDF4 file (there should always be a match)
        slc[self._axes['time']] = match(xr_time, self._nc_time)
        return self._nc_data.variables[self._varname][slc]

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
