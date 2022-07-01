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
        """Model input data.

        Parameters
        ----------
        model: hm.Model
            Model object which inherits from `hm.Model`.
        filename: str
            Filename of the dataset.
        nc_varname: str 
            Variable name in netCDF file.
        model_varname: str
            Model variable name.
        is_1d: bool, optional
            Whether the dataset has a one-dimensional 
            representation of space (i.e. a set of points). This 
            is opposed to a two-dimensional representation which will 
            have coordinates to identify the location of each point.
        xy_dimname: str, optional
            If the dataset is one-dimensional, this parameter specifies the name
            of the space dimension, which is often non-standard (e.g. 'land').
        factor: float, optional
            Factor by which to multiply data values.
        offset: float, optional
            Offset to add to data values.
        """
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

    def update(self, method, **kwargs):
        vars(self.model)[self.dataset_varname].select(
            time=self.model.time.curr_time, method=method, **kwargs
        )
        vars(self.model)[self.model_varname] = (vars(self.model)[self.dataset_varname].values * self.factor) + self.offset

    def dynamic(self, method='nearest'):
        # TODO: api to HmInputData to work out whether temporal or not
        if vars(self.model)[self.dataset_varname].is_temporal:
            self.update(method)


class HmSpaceInputData(HmInputData):
    def dynamic(self):
        pass
