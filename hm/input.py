#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np

from .api import open_hmdataarray

# Wishlist:
# * update at intervals greater than the model timestep (e.g. annually)

class HmInputData(object):
    def __init__(
            self,
            model,
            filename,
            nc_varname,
            model_varname
    ):
        self.model = model
        self.filename = filename
        self.nc_varname = nc_varname
        self.model_varname = model_varname        

    def initial(self):
        self.read()
        
    def read(self):
        vars(self.model)[self.model_varname] = open_hmdataarray(
            self.filename,
            self.nc_varname,
            self.model.domain
        )
        
    def update(self, method):
        vars(self.model)[self.model_varname].select(
            time=self.model.time.curr_time, method=method
        )
        # vars(self.model)[self.model_varname].select(
        #     time=self.model.time.curr_time, method='nearest'
        # )
        
    def dynamic(self, method='nearest'):
        self.update(method)
        
class HmSpaceInputData(object):
    def dynamic(self):
        pass
