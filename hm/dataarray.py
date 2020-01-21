#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import math
import string
import datetime
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
        self._dims = get_dimension_names(dataarray_or_dataset, is_1d, xy_dimname)
        self._coords = get_coordinates(dataarray_or_dataset, self._dims)
        self._is_1d = is_1d & ('xy' in self.dims)
        self._is_2d = (not self._is_1d) & all(dim in self.dims for dim in ('x','y'))
        self._is_spatial = self._is_1d | self._is_2d
        self.spatial_extent()
        
    def spatial_extent(self):
        if self.is_2d:
            self._extent = spatial_extent(self._coords)
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
        super().__init__(dataarray, is_1d, xy_dimname)
        self._domain = domain
        self.subset()
        
    def subset(self):
        self.get_index()
        self._data = self._data.sel(self.index)
        
    def get_index(self):
        index_dict = {}
        for dim in self.dims:
            index_dict[dim] = self._domain.coords[dim]        
        self.index = index_dict
    def select(self):
        pass    
    def load(self):
        self._data.load()
    
class HmSpaceTimeDataArray(HmSpaceDataArray):
    def get_index(self):
        super().get_index()
        time_slc = slice(self._domain.starttime, self._domain.endtime)
        self.index['time'] = time_slc        
    def select(self, time):
        # 'time' may be a slice, as in slice(starttime, endtime)
        return self._data.sel({self._dims['time'] : time})
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
