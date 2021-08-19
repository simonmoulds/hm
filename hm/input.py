#!/usr/bin/env python
# -*- coding: utf-8 -*-

from .api import open_hmdataarray

# Wishlist:
# * update at intervals greater than the model timestep (e.g. annually)


class HmInputData(object):
    def __init__(
            self,
            model,
            filename,
            nc_varname,
            model_varname,
            is_1d=False,
            xy_dimname=None,
            factor=1.,
            offset=0.
    ):
        self.model = model
        self.filename = filename
        self.nc_varname = nc_varname
        self.dataset_varname = model_varname + '_dataset'
        self.model_varname = model_varname
        self.is_1d = is_1d
        self.xy_dimname = xy_dimname
        self.factor = float(factor)
        self.offset = float(offset)

    def initial(self):
        self.read()

    def read(self):
        # TODO: think about whether it's really necessary for model
        # developers to be able to access the underlying dataset,
        # or whether they should simply interface with it through this class
        vars(self.model)[self.dataset_varname] = open_hmdataarray(
            self.filename,
            self.nc_varname,
            self.model.domain,
            self.is_1d,
            self.xy_dimname,
            self.model.is_1d
        )

    def update(self, method):
        vars(self.model)[self.dataset_varname].select(
            time=self.model.time.curr_time, method=method
        )
        vars(self.model)[self.model_varname] = (vars(self.model)[self.dataset_varname].values * self.factor) + self.offset

    def dynamic(self, method='nearest'):
        # TODO: api to HmInputData to work out whether temporal or not
        if vars(self.model)[self.dataset_varname].is_temporal:
            self.update(method)


class HmSpaceInputData(HmInputData):
    def dynamic(self):
        pass
