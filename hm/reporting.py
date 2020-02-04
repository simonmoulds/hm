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
    all_vars = sorted(set(all_vars))
    var_dict['all'] = all_vars
    return box.Box(var_dict, frozen_box=True)

def _get_shortname(varname):
    try:
        shortname = variable_list.netcdf_shortname[varname]
    except KeyError:
        pass
    return shortname

def _get_standard_name(varname):
    try:
        standard_name = variable_list.netcdf_standard_name[varname]
    except KeyError:
        standard_name = variable_list.netcdf_shortname[varname]
    return standard_name

def _get_long_name(varname):
    try:
        long_name = variable_list.netcdf_long_name[varname]
    except KeyError:
        long_name = variable_list.netcdf_shortname[varname]
    return long_name

def _get_dimensions(varname):
    try:
        dimensions = variable_list.netcdf_dimensions[varname]
    except KeyError:
        pass
    return dimensions

def _get_units(varname):
    try:
        units = variable_list.netcdf_units[varname]
    except:
        pass
    return units

def _get_description(varname):
    try:
        description = variable_list.netcdf_description[varname]
    except KeyError:
        description = None
    return description    

def _get_calendar(varname):
    try:
        calendar = variable_list.netcdf_calendar[varname]
    except KeyError:
        calendar = None
    return calendar

def _get_variable_datatype(varname):
    try:
        datatype = variable_list.netcdf_datatype[varname]
    except:
        datatype = 'f4'
    return datatype

def _get_variable_attributes(varname):
    attr = {
        'shortname'     : _get_shortname(varname),
        'standard_name' : _get_standard_name(varname),
        'long_name'     : _get_long_name(varname),
        'units'         : _get_units(varname),
        'calendar'      : _get_calendar(varname),
        'description'   : _get_description(varname),
        'dimensions'    : _get_dimensions(varname),
        'datatype'      : _get_datatype(varname)
    }
    return box.Box(attr, frozen_box=True)

class _netcdf(object):
    def __init__(self, varname, filename):
        self.varname = varname
        self.filename = filename
        self.attr = _get_variable_attributes(self.varname)
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
        for dimname in self.dimensions:            
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
        attr = _get_variable_attributes(dimname)
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
        var.calendar = attr.units
        self.netcdf_obj.sync()

    def add_nontime_dimension(self, dimname, **kwargs):
        """Add non-time dimension to a netCDF file."""
        attr = _get_variable_attributes(dimname)
        self.netcdf_obj.createDimension(attr.shortname, None)
        if dim in allowed_x_dim_names + allowed_y_dim_names:
            # if spatial dimension, ensure sufficient precision
            # to deal with resolutions with repeating decimals
            # (e.g. 1/12 degree)
            kwargs['zlib'] = True
            kwargs['least_significant_digit'] = 16            
        var = self.netcdf_obj.createVariable(
            shortname,
            datatype,
            dimensions,
            **kwargs
        )
        var.standard_name = variable_list.netcdf_standard_name[dimname]
        var.long_name = variable_list.netcdf_long_name[dimname]
        var.units = variable_list.netcdf_unit[dimname]
        var[:] = np.array(self._model.coords[dimname])            
        self.netcdf_obj.sync()

    def get_time_axis(self):
        if self.is_temporal:
            self.time_axis = [index for index, value in enumerate(self.dimensions) if value in allowed_t_dim_names][0]
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
            self.add_time_varying_data()
        else:
            self.add_static_data()
        self.netcdf_obj.sync()

    def add_time_varying_data(self, data):
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
        slc = [slice(None)] * len(self.dimensions)
        slc[self.time_axis] = time_index
        self.netcdf_obj.variables[self.attr.shortname][slc] = data
        
    def add_static_data(self, data):
        self.netcdf_obj.variables[self.attr.shortname][:] = data
        
class SummaryVariable(object):
    def __init__(
            self,
            model,
            varname,
            freq
    ):        
        self._model = model
        self.varname = varname
        self.filename = self.make_filename()
        self._netcdf = _netcdf(self.varname, self.filename)
        self.data = np.zeros_like(
            vars(model)[varname].shape[:-1] + (self._model.domain.mask.shape)
        )
        self.get_reporting_times(freq=freq)
        self.n_timestep = 0

    def get_reporting_times(self, **kwargs):
        # TODO: use consistent time classes - pandas, datetime, numpy etc.
        # I think pandas is the best option - seems to focus on time series, which fits nicely with hydrological modelling 
        # https://pandas.pydata.org/pandas-docs/stable/user_guide/timeseries.html
        self.reporting_times = pd.date_range(
            self._model.model_time.starttime,
            self._model.model_time.endtime,
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
        self.filename = os.path.join(
            self._model.config.output_directory,
            self._model.config.output_prefix + '_output_' + self.varname + '.nc4'
        )

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
        self.data += vars(self._model)[self.varname][self._model.domain.mask]
        self.n_timestep += 1
        if self._model.model_time.curr_time == self._end_of_current_reporting_period:
            self.data /= self.n_timestep
            self.to_netcdf()
            self.reset()

class MaxSummaryVariable(SummaryVariable):
    def update(self):
        self.data = self.data.clip(vars(self._model)[self.varname], None)
        self.n_timestep += 1
        if self._model.model_time.curr_time == self._end_of_current_reporting_period:
            self.to_netcdf()
            self.reset()
            
class MinSummaryVariable(SummaryVariable):
    def update(self):
        self.data = self.data.clip(None, vars(self._model)[self.varname])
        self.n_timestep += 1
        if self._model.model_time.curr_time == self._end_of_current_reporting_period:
            self.to_netcdf()
            self.reset()

class EndSummaryVariable(SummaryVariable):
    def update(self):
        self.n_timestep += 1
        if self._model.model_time.curr_time == self._end_of_current_reporting_period:
            self.data = vars(self._model)[self.varname]
            self.to_netcdf()
            self.reset()

class TotalSummaryVariable(SummaryVariable):
    def update(self):
        self.data += vars(self._model)[self.varname]
        self.n_timestep += 1
        if self._model.model_time.curr_time == self._end_of_current_reporting_period:
            self.to_netcdf()
            self.reset()

class DummyReporting(object):
    def __init__(self):
        pass

    def report(self):
        pass
            
class Reporting(object):
    def __init__(
            self,
            model,
            reporting_options=None,  # not used currently
            run_id=None              # not used currently
    ):
        self._model = model
        self._outdir = self._model._config.output_dir
        self._reporting_variables = _get_reporting_variables(self._model.config)
        self.initiate_reporting()

    def create_summary_variables(self):
        output_dict = {}
        for option, varnames in self._reporting_variables.items():
            if varnames is not None:
                for var in varnames:
                    output_dict[var + '_' + option] = SummaryVariable(self._model, var, freq)

    def update(self):
        for key, value in output_dict.items():
            value.update()
            
