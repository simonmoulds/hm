#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import shutil
import numpy as np
import datetime
import box

from .constants import *
from .utils import *

import logging
logger = logging.getLogger(__name__)

def get_variable_names_for_reporting(config, option):
    try:
        return [str(var.strip()) for var in config.REPORTING[option].split(',')]
    except KeyError:
        return []

def _get_reporting_variables(config):
    # var_dict = {'all' : []}
    var_dict = {}
    all_vars = []
    for option in allowed_reporting_options:
        var_names = get_variable_names_for_reporting(config, option)
        var_dict[option] = var_names
        all_vars += var_names
    # all_vars = sorted(set(all_vars))
    # var_dict['all'] = all_vars
    return box.Box(var_dict, frozen_box=True)

def _get_shortname(varname, variable_list):
    try:
        shortname = variable_list.netcdf_shortname[varname]
    except KeyError:
        pass
    return shortname

def _get_standard_name(varname, variable_list):
    try:
        standard_name = variable_list.netcdf_standard_name[varname]
    except KeyError:
        standard_name = variable_list.netcdf_shortname[varname]
    return standard_name

def _get_long_name(varname, variable_list):
    try:
        long_name = variable_list.netcdf_long_name[varname]
    except KeyError:
        long_name = variable_list.netcdf_shortname[varname]
    return long_name

def _get_dimensions(varname, variable_list):
    try:
        dimensions = variable_list.netcdf_dimensions[varname]
    except KeyError:
        pass
    return dimensions

def _get_units(varname, variable_list):
    try:
        units = variable_list.netcdf_units[varname]
    except:
        pass
    return units

def _get_description(varname, variable_list):
    try:
        description = variable_list.netcdf_description[varname]
    except KeyError:
        description = None
    return description    

def _get_calendar(varname, variable_list):
    try:
        calendar = variable_list.netcdf_calendar[varname]
    except KeyError:
        calendar = None
    return calendar

def _get_datatype(varname, variable_list):
    try:
        datatype = variable_list.netcdf_datatype[varname]
    except:
        datatype = 'f4'
    return datatype

def _get_variable_attributes(varname, variable_list):
    attr = {
        'shortname'     : _get_shortname(varname, variable_list),
        'standard_name' : _get_standard_name(varname, variable_list),
        'long_name'     : _get_long_name(varname, variable_list),
        'units'         : _get_units(varname, variable_list),
        'calendar'      : _get_calendar(varname, variable_list),
        'description'   : _get_description(varname, variable_list),
        'dimensions'    : _get_dimensions(varname, variable_list),
        'datatype'      : _get_datatype(varname, variable_list)
    }
    return box.Box(attr, frozen_box=True)

class _netcdf(object):
    def __init__(self, model, varname, filename, variable_list):
        self.model = model
        self.varname = varname
        self.filename = filename
        self.variable_list = variable_list
        self.attr = _get_variable_attributes(self.varname, self.variable_list)
        self.netcdf_obj = open_netcdf(self.filename, mode='w')
        self.add_global_attributes()    
        self.add_dimensions()
        self.get_time_axis()
        self.add_variable()

    def add_global_attributes(self):
        """Add global netCDF attributes."""
        try:
            for name, value in global_netcdf_attributes.items():
                self.netcdf_obj.setncattr(name, value)
        except:
            pass
        
    def add_dimensions(self):
        """Add dimensions to netCDF object."""
        # TODO: what is the difference between a dimension and coordinate
        # See discussion: https://math.stackexchange.com/questions/3327858/terminology-dimension-vs-coordinate
        is_temporal = False
        # for dimname in self.dimensions:
        for dimname in self.attr.dimensions:
            if dimname in allowed_t_dim_names:
                # i.e. time dimension, which should be unlimited
                self.add_time_dimension(dimname)
                is_temporal = True
            else:
                self.add_nontime_dimension(dimname)
        self.is_temporal = is_temporal
        
    def add_time_dimension(self, dimname, **kwargs):
        """Add time dimension to a netCDF file.

        Time dimension is unlimited so that data can 
        be added during simulation.
        """
        # get dimension attributes, then create dimension
        attr = _get_variable_attributes(dimname, self.variable_list)
        self.netcdf_obj.createDimension(attr.shortname, None)
        var = self.netcdf_obj.createVariable(
            attr.shortname,
            attr.datatype,
            attr.dimensions,
            **kwargs
        )
        # add additional variable attributes
        var.standard_name = attr.standard_name
        var.long_name = attr.long_name
        var.units = attr.units
        var.calendar = attr.calendar
        self.netcdf_obj.sync()

    def add_nontime_dimension(self, dimname, **kwargs):
        """Add non-time dimension to a netCDF file."""
        attr = _get_variable_attributes(dimname, self.variable_list)
        self.netcdf_obj.createDimension(attr.shortname, None)
        if dimname in allowed_x_dim_names + allowed_y_dim_names:
            # if spatial dimension, ensure sufficient precision
            # to deal with resolutions with repeating decimals
            # (e.g. 1/12 degree)
            kwargs['zlib'] = True
            kwargs['least_significant_digit'] = 16            
        var = self.netcdf_obj.createVariable(
            attr.shortname,
            attr.datatype,
            attr.dimensions,
            **kwargs
        )
        var.standard_name = self.variable_list.netcdf_standard_name[dimname]
        var.long_name = self.variable_list.netcdf_long_name[dimname]
        var.units = self.variable_list.netcdf_units[dimname]
        if dimname in allowed_x_dim_names:
            hm_dimname = 'x'
        elif dimname in allowed_y_dim_names:
            hm_dimname = 'y'
        elif dimname in allowed_z_dimnames:
            hm_dimname = 'z'
        var[:] = np.array(self.model.domain.coords[hm_dimname])
        self.netcdf_obj.sync()

    def get_time_axis(self):
        if self.is_temporal:
            self.time_axis = [index for index, value in enumerate(self.model.domain._dims) if value in allowed_t_dim_names][0]
        else:
            self.time_axis = None

    def add_variable(self, **kwargs):
        var = self.netcdf_obj.createVariable(
            self.attr.shortname,
            self.attr.datatype,
            self.attr.dimensions,
            **kwargs
        )
        var.standard_name = self.attr.standard_name
        var.long_name = self.attr.long_name
        var.units = self.attr.units
        self.netcdf_obj.sync()
        
    def add_data(self, data):
        """Add data to the netCDF file."""
        if self.is_temporal:
            self.add_time_varying_data(data)
        else:
            self.add_static_data(data)
        self.netcdf_obj.sync()

    def add_time_varying_data(self, data):
        print(data.shape)
        timevar = self.netcdf_obj.variables['time']
        try:
            time_index = len(timevar)
        except:
            time_index = 0
        timevar[time_index] = nc.date2num(
            self.model.time.timestamp,
            timevar.units,
            timevar.calendar
        )        
        slc = [slice(None)] * len(self.attr.dimensions)
        slc[self.time_axis] = time_index
        self.netcdf_obj.variables[self.attr.shortname][slc] = data
        
    def add_static_data(self, data):
        self.netcdf_obj.variables[self.attr.shortname][:] = data
        
class SummaryVariable(object):
    def __init__(
            self,
            model,
            varname,
            variable_list,
            freq
    ):        
        self.model = model
        self.varname = varname
        # self.filename = self.make_filename()
        self.filename = 'test.nc'
        # print(self.filename)
        self._netcdf = _netcdf(
            self.model,
            self.varname,
            self.filename,
            variable_list
        )
        # self.data = np.zeros_like(
        #     vars(model)[varname].shape[:-1] + (self.model.domain.mask.shape)
        # )
        self.get_reporting_times(freq=freq)
        self._end_of_current_reporting_period = self.reporting_times[0]
        self._counter = 0
        self.n_timestep = 0
        
    def get_reporting_times(self, **kwargs):
        # TODO: use consistent time classes - pandas, datetime, numpy etc.
        # I think pandas is the best option - seems to focus on time series, which fits nicely with hydrological modelling 
        # https://pandas.pydata.org/pandas-docs/stable/user_guide/timeseries.html
        self.reporting_times = pd.date_range(
            self.model.time.starttime,
            self.model.time.endtime,
            **kwargs
        )# .to_pydatetime().tolist()
        # if reporting_times[-1] < self._model.model_time.endtime:
        #     reporting_times[-1] = self._model.model_time.endtime
        # TODO: data type conversion
        # self.reporting_times = reporting_times        
        if len(self.reporting_times) == 0:
            warnings.warn(
                'Summary period ' + freq + ' is longer than '
                'the model run duration: not creating output '
                'for variable ' + varname + ' at frequency '
                + freq
            )
        
    def make_filename(self):
        """Make the filename."""
        # TODO: allow user to specify preferred netCDF extension
        self.filename = 'test.nc'# os.path.join(
        #     self.model.config.output_directory,
        #     self.model.config.output_prefix + '_output_' + self.varname + '.nc4'
        # )
    # def initial(self):
    #     self.data = np.zeros_like(
    #         vars(self.model)[self.varname].shape[:-1] + (self.model.domain.mask.shape)
    #     )
    #     print(self.data.shape)        
    def update(self):
        pass
            
    def to_netcdf(self):
        self._netcdf.add_data(self.data)

    def reset(self):
        self.data.fill(0)
        self.n_timestep = 0
        self._counter += 1
        self._end_of_current_reporting_period = self.reporting_times[self._counter]

class MeanSummaryVariable(SummaryVariable):
    def update(self):
        # if self.model.time.is_first_timestep:
        #     self.data = np.zeros_like(
        #         vars(self.model)[self.varname].shape[:-1] + (self.model.domain.mask.shape)
        #     )            
        self.data += vars(self.model)[self.varname][self.model.domain.mask]
        self.n_timestep += 1
        if self.model.time.curr_time == self._end_of_current_reporting_period:
            self.data /= self.n_timestep
            self.to_netcdf()
            self.reset()

class MaxSummaryVariable(SummaryVariable):
    def update(self):
        self.data = self.data.clip(vars(self.model)[self.varname], None)
        self.n_timestep += 1
        if self.model.time.curr_time == self._end_of_current_reporting_period:
            self.to_netcdf()
            self.reset()
            
class MinSummaryVariable(SummaryVariable):
    def update(self):
        self.data = self.data.clip(None, vars(self.model)[self.varname])
        self.n_timestep += 1
        if self.model.time.curr_time == self._end_of_current_reporting_period:
            self.to_netcdf()
            self.reset()

class EndSummaryVariable(SummaryVariable):
    def update(self):
        self.n_timestep += 1
        if self.model.time.curr_time == self._end_of_current_reporting_period:
            self.data = vars(self.model)[self.varname]
            self.to_netcdf()
            self.reset()

class TotalSummaryVariable(SummaryVariable):
    def update(self):
        if self.model.time.is_first_timestep:
            self.data = np.zeros(
                vars(self.model)[self.varname].shape[:-1] + (self.model.domain.mask.shape),
                dtype=np.float64
            )
        self.data[self.model.domain.mask] += vars(self.model)[self.varname]
        self.n_timestep += 1
        if self.model.time.curr_time == self._end_of_current_reporting_period:
            self.to_netcdf()
            self.reset()

class DummyReporting(object):
    def __init__(self):
        pass

    def update(self):
        pass

def open_summary_variable(model, varname, option, variable_list):
    if option == 'daily_total':
        return TotalSummaryVariable(model, varname, variable_list, freq='1D')

class Reporting(object):
    def __init__(
            self,
            model,
            variable_list,
            reporting_options=None,  # not used currently
            run_id=None              # not used currently
    ):
        self.model = model
        self.variable_list = variable_list
        self.outdir = self.model.config.output_directory
        self.reporting_variables = _get_reporting_variables(self.model.config)
        self.create_summary_variables()

    def create_summary_variables(self):
        self.output_variables = {}
        for option, varnames in self.reporting_variables.items():
            if varnames is not None:
                for varname in varnames:
                    # print(option, varname)
                    print(option,varname)
                    self.output_variables[varname + '_' + option] = open_summary_variable(self.model, varname, option, self.variable_list)
                    # self.output_variables[var + '_' + option] = SummaryVariable(self.model, var, freq)

    def initial(self):
        pass
    
    def update(self):
        for _, value in self.output_variables.items():
            value.update()
