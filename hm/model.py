#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import shutil
import sys
import math
import gc
import numpy as np 
import rasterio
import xarray as xr

class Model(object):
    def __init__(
            self,
            config,
            time,
            init=None
    ):
        self._config = config
        self._time = modeltime
        self.set_domain()
    def set_domain(self):
        self._domain = set_domain(self._config.mask, self._time)
        
# class Model(object):
    
#     def __init__(self, configuration, modelTime, initialState = None):
#         self._configuration = configuration
#         self._modelTime = modelTime
#         self.set_clone_map()
#         self.set_landmask()
#         self.set_grid_cell_area()
#         self.get_model_dimensions()
        
#     def set_clone_map(self):
#         # different between 1D/2D: 1D use netCDF; 2D use netCDF or GDAL
#         self.cloneMapFileName = str(self._configuration.cloneMap)
#         # self.cloneMap = file_handling.read_GDAL(self.cloneMapFileName)
#         try:
#             cloneVarName = 'clone'
#             ds = xr.open_dataset(self.cloneMapFileName)
#             self.clone = ds.clone.values
#             self.clone_x_coord = ds.lon.values
#             self.clone_y_coord = ds.lat.values
#             ds.close()
#         except:
#             ds = xr.open_rasterio(self.cloneMapFileName)
#             self.clone = ds.values
#             self.clone_x_coord = ds.x.values
#             self.clone_y_coord = ds.y.values
#             ds.close()
            
#     def set_landmask(self):
#         # different between 1D/2D: 1D use netCDF; 2D use netCDF or GDAL
#         try:
#             self.landmask = file_handling.gdal_to_array(self._configuration.landmask, cloneMapFileName=self.cloneMapFileName)
#         except:
#             self.landmask = self.clone.copy()
#             self.landmask[:] = 1
#             # self.landmask = file_handling.read_map_clone(self.configuration.landmask, self.cloneMap)
#         self.landmask = self.landmask > 0

#     def set_grid_cell_area(self):
#         # only Model2D: meaningless for point runs
#         # TEMPORARY: need to make this more robust!!!
#         try:
#             grid_cell_area = file_handling.netcdf_to_arrayWithoutTime(
#                 str(self._configuration.MODEL_GRID['gridCellAreaInputFile']),
#                 str(self._configuration.MODEL_GRID['gridCellAreaVariableName']),
#                 cloneMapFileName = self.cloneMapFileName
#             )
#             grid_cell_area = np.float64(grid_cell_area)
#             self.grid_cell_area = grid_cell_area[self.landmask]
            
#         except:
            
#             try:
#                 grid_cell_area = file_handling.gdal_to_array(
#                     str(self._configuration.MODEL_GRID['gridCellAreaInputFile']),
#                     cloneMapFileName=self.cloneMapFileName
#                     )
#                 self.grid_cell_area = grid_cell_area[self.landmask]
                
#             except:                
#                 self.grid_cell_area = np.ones_like(self.landmask)[self.landmask]

#     def get_model_dimensions(self):
#         # different between 1D/2D
#         """Function to set model dimensions"""
#         self.latitudes = self.clone_y_coord
#         self.nLat = self.clone_y_coord.size
#         # self.latitudes = self.cloneMap.y.values
#         # self.nLat = self.cloneMap.y.size
#         self.longitudes = self.clone_x_coord
#         self.nLon = self.clone_x_coord.size
#         # self.longitudes = self.cloneMap.x.values
#         # self.nLon = self.cloneMap.x.size
#         self.nCell = int(np.sum(self.landmask))
#         self.dimensions = {
#             'time'     : None,
#             'lat'      : self.latitudes,
#             'lon'      : self.longitudes,
#         }
        
#     @property
#     def configuration(self):
#         return self._configuration
