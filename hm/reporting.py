#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import numpy as np
import box
import xarray

from .constants import *
from .utils import *

import logging
logger = logging.getLogger(__name__)


def get_variable_names_for_reporting(config, section, option):
    try:
        varnames = vars(config)[section][option]
        if type(varnames) is list:
            return varnames
        else:
            return [varnames]
    except KeyError:
        return []


def _get_summary_variables(config,  section):
    var_dict = {}
    for option in allowed_reporting_options:
        sum_var_names = get_variable_names_for_reporting(
            config, section, option + '_summary'
        )
        if len(sum_var_names) > 0:
            var_dict[option] = sum_var_names
    return box.Box(var_dict, frozen_box=True)

def _get_reporting_variables(config, section):
    var_dict = {}
    for option in allowed_reporting_options:
        var_names = get_variable_names_for_reporting(
            config, section, option
        )
        var_dict[option] = var_names
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
    """Retrieve netCDF attributes for model variable."""
    attr = {
        'shortname': _get_shortname(varname, variable_list),
        'standard_name': _get_standard_name(varname, variable_list),
        'long_name': _get_long_name(varname, variable_list),
        'units': _get_units(varname, variable_list),
        'calendar': _get_calendar(varname, variable_list),
        'description': _get_description(varname, variable_list),
        'dimensions': _get_dimensions(varname, variable_list),
        'datatype': _get_datatype(varname, variable_list)
    }
    return box.Box(attr, frozen_box=False)


def _get_standard_dimname(dimname):
    if dimname in allowed_xy_dim_names:
        return 'xy'
    if dimname in allowed_x_dim_names:
        return 'x'
    elif dimname in allowed_y_dim_names:
        return 'y'
    elif dimname in allowed_z_dim_names:
        return 'z'
    else:
        return dimname


class _netcdf(object):
    
    def __init__(
            self,
            model,
            varname,
            filename,
            variable_list  # ,
            # transpose = False
    ):
        self.model = model
        self.varname = varname
        self.filename = filename
        self.ncdf = open_netcdf(self.filename, mode='w')
        self.variable_list = variable_list
        # self.transpose = bool(transpose)
        self.attr = _get_variable_attributes(self.varname, self.variable_list)
        if self.model.domain.is_1d:
            self.attr.dimensions = self.attr.dimensions + ('space',)
        else:
            self.attr.dimensions = self.attr.dimensions + ('lat', 'lon')
        # ensure attr.dimensions are unique
        self.attr.dimensions = list(dict.fromkeys(self.attr.dimensions))            
        self.add_global_attributes()
        self.add_dimensions()
        self.get_time_axis()
        self.add_variable(self.attr)

    def add_global_attributes(self):
        try:
            for name, value in global_netcdf_attributes.items():
                self.ncdf.setncattr(name, value)
        except:
            pass

    def add_dimensions(self):
        """Add dimensions to netCDF object."""
        # TODO: what is the difference between a dimension and coordinate
        # See discussion: https://math.stackexchange.com/questions/3327858/terminology-dimension-vs-coordinate

        # if self.model.domain.is_1d:
        #     self.attr.dimensions = self.attr.dimensions + ('space',)
        # else:
        #     self.attr.dimensions = self.attr.dimensions + ('lat', 'lon')
        # # ensure attr.dimensions are unique
        # self.attr.dimensions = list(dict.fromkeys(self.attr.dimensions))        
        is_temporal = False
        for dim in self.attr.dimensions:
            if dim in allowed_t_dim_names:
                self.add_time_dimension(dim)
                is_temporal = True
            else:
                if self.model.domain.is_1d:
                    if not (dim in allowed_x_dim_names) | (dim in allowed_y_dim_names):
                        self.add_nontime_dimension(dim)
                else:
                    self.add_nontime_dimension(dim)

        # TODO: this doesn't work for 1D outputs, because the dimensions in
        # variable_list include lat/lon
        # - Remove lat/lon from variable_list dimensions
        # - Test if 1D or 2D
        # - Test if lat/long
        # - If 1D, we need to write spatial coordinates as variables
        if self.model.domain.is_1d:
            for dimname in ['lat','lon']:
                attr = _get_variable_attributes(dimname, self.variable_list)
                attr.dimensions = ('space',)
                self.add_variable(attr)
                standard_dimname = _get_standard_dimname(dimname)
                self.ncdf.variables[dimname][:] = np.array(self.model.domain.coords[standard_dimname])
                self.ncdf.sync()
            
        self.is_temporal = is_temporal

    def add_time_dimension(self, dimname, **kwargs):
        """Add time dimension to a netCDF file.

        Time dimension is unlimited so that data can
        be added during simulation.
        """
        attr = _get_variable_attributes(dimname, self.variable_list)
        self.ncdf.createDimension(attr.shortname, None)
        var = self.ncdf.createVariable(
            attr.shortname,
            attr.datatype,
            attr.dimensions,
            **kwargs
        )
        var.standard_name = attr.standard_name
        var.long_name = attr.long_name
        var.units = attr.units
        var.calendar = attr.calendar
        self.ncdf.sync()

    def add_nontime_dimension(self, dimname, **kwargs):
        """Add non-time dimension to a netCDF file."""
        attr = _get_variable_attributes(dimname, self.variable_list)
        # since nontime dimensions are limited (unlike time, which
        # is conventionally unlimited), we need to get the dimension
        # size
        standard_dimname = _get_standard_dimname(dimname)
        dimvals = np.array(self.model.domain.coords[standard_dimname])
        self.ncdf.createDimension(attr.shortname, len(dimvals))
        if dimname in allowed_x_dim_names + allowed_y_dim_names:
            # if spatial dimension, ensure sufficient precision
            # to deal with resolutions with repeating decimals
            # (e.g. 1/12 degree)
            kwargs['zlib'] = True
            kwargs['least_significant_digit'] = 16
            
        var = self.ncdf.createVariable(
            attr.shortname,
            attr.datatype,
            attr.dimensions,
            **kwargs
        )
        var.standard_name = self.variable_list.netcdf_standard_name[dimname]
        var.long_name = self.variable_list.netcdf_long_name[dimname]
        var.units = self.variable_list.netcdf_units[dimname]
        standard_dimname = _get_standard_dimname(dimname)
        var[:] = dimvals
        self.ncdf.sync()

    def get_time_axis(self):
        if self.is_temporal:
            self.time_axis = [index for index, value in enumerate(
                self.attr.dimensions) if value in allowed_t_dim_names][0]
            # self.time_axis = [index for index, value in enumerate(self.model.domain._dims) if value in allowed_t_dim_names][0]
        else:
            self.time_axis = None

    def add_variable(self, attr, **kwargs):
        """Add variable to netCDF file."""
        var = self.ncdf.createVariable(
            attr.shortname,
            attr.datatype,
            attr.dimensions,
            **kwargs
        )
        try:
            var.standard_name = attr.standard_name
        except:
            pass
        try:
            var.long_name = attr.long_name
        except:
            pass
        try:
            var.units = attr.units
        except:
            pass
        self.ncdf.sync()

    def add_data(self, data):
        """Add data to netCDF file."""
        if self.is_temporal:
            self.add_time_varying_data(data)
        else:
            self.add_static_data(data)
        self.ncdf.sync()

    def add_time_varying_data(self, data):
        """Add time-varying data to netCDF file."""
        timevar = self.ncdf.variables['time']
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
        self.ncdf.variables[self.attr.shortname][slc] = data

    def add_static_data(self, data):
        """Add time-invariant data to the netCDF file."""
        self.ncdf.variables[self.attr.shortname][:] = data


class ReportingVariable(object):
    def __init__(
            self,
            model,
            varname,
            sample,
            variable_list,
            freq,
            suffix  # ,
            # transpose = False
    ):
        self.model = model
        self.varname = varname
        self.filename = os.path.join(
            self.model.config.output_directory,
            'netcdf',
            'hm_output_' + suffix + '_' + self.varname + '.' + str(sample).zfill(3) + '.nc'
        )        
        # self.make_filename(suffix, sample)
        # self.transpose = bool(transpose)
        self._netcdf = _netcdf(
            self.model,
            self.varname,
            self.filename,
            variable_list  # ,
            # self.transpose
        )
        self.get_reporting_times(freq=freq)
        self.n_reporting_times = len(self.reporting_times)
        try:
            self.end_of_current_reporting_period = self.reporting_times[0]
        except:
            self.end_of_current_reporting_period = self.model.time.endtime
        self.counter = 0
        self.n_timestep = 0

    def initial(self):
        var = vars(self.model)[self.varname]
        # if self.transpose:
        #     var = var.T
        if self.model.is_1d:
            shape = var.shape[:-1]
        else:
            shape = var.shape[:-2]
        self.data = np.zeros(
            shape + (self.model.domain.mask.shape), dtype=np.float64)
        
    def get_reporting_times(self, **kwargs):
        self.reporting_times = pd.date_range(
            self.model.time.starttime,
            self.model.time.endtime,
            **kwargs
        )
        self.reporting_times = self.reporting_times.unique()
        if len(self.reporting_times) == 0:
            warnings.warn(
                'Reporting period ' + freq + ' is longer than '
                'the model run duration: not creating output '
                'for variable ' + varname + ' at frequency '
                + freq
            )

    # def make_filename(self, suffix, sample):
    #     self.filename = os.path.join(
    #         self.model.config.output_directory,
    #         'netcdf',
    #         'hm_output_' + suffix + '_' + self.varname + '.' + str(sample).zfill(3) + '.nc'
    #     )

    def to_netcdf(self):
        self._netcdf.add_data(self.data)

    def reset(self):
        self.data.fill(0)
        self.n_timestep = 0
        self.counter += 1
        if self.counter < self.n_reporting_times:
            self.end_of_current_reporting_period = self.reporting_times[self.counter]
        else:
            self.end_of_current_reporting_period = self.model.time.endtime
            
    def update(self):
        pass

class MeanReportingVariable(ReportingVariable):
    def update(self):
        self.data[...,
                  self.model.domain.mask] += vars(self.model)[self.varname]
        self.n_timestep += 1
        if self.model.time.curr_time == self.end_of_current_reporting_period:
            self.data /= self.n_timestep
            self.to_netcdf()
            self.reset()


class MaxReportingVariable(ReportingVariable):
    def update(self):
        data = self.data[..., self.model.domain.mask]
        data = data.clip(vars(self.model)[self.varname], None)
        self.data[..., self.model.domain.mask] = data
        self.n_timestep += 1
        if self.model.time.curr_time == self.end_of_current_reporting_period:
            self.to_netcdf()
            self.reset()


class MinReportingVariable(ReportingVariable):
    def update(self):
        data = self.data[..., self.model.domain.mask]
        data = data.clip(None, vars(self.model)[self.varname])
        self.dataself.data[..., self.model.domain.mask] = data
        self.n_timestep += 1
        if self.model.time.curr_time == self.end_of_current_reporting_period:
            self.to_netcdf()
            self.reset()


class EndReportingVariable(ReportingVariable):
    def update(self):
        self.n_timestep += 1
        if self.model.time.curr_time == self.end_of_current_reporting_period:
            self.data[..., self.model.domain.mask] = vars(self.model)[
                self.varname]
            self.to_netcdf()
            self.reset()


class TotalReportingVariable(ReportingVariable):
    def update(self):
        self.data[...,
                  self.model.domain.mask] += vars(self.model)[self.varname]
        self.n_timestep += 1
        if self.model.time.curr_time == self.end_of_current_reporting_period:
            self.to_netcdf()
            self.reset()


class DummyReporting(object):
    def __init__(self):
        pass

    def initial(self):
        pass

    def dynamic(self):
        pass


def open_reporting_variable(model, varname, option, sample, variable_list):
    # TODO - complete for all frequencies/statistics
    if option == 'daily_total':
        return TotalReportingVariable(model, varname, sample, variable_list, freq='1D', suffix='daily_total')
    if option == 'year_max':
        return MaxReportingVariable(model, varname, sample, variable_list, freq='1Y', suffix='year_max')

# OLD (no sample):

# def open_reporting_variable(model, varname, option, variable_list):
#     # TODO - complete for all frequencies/statistics
#     if option == 'daily_total':
#         return TotalReportingVariable(model, varname, variable_list, freq='1D', suffix='daily_total')
#     if option == 'year_max':
#         return MaxReportingVariable(model, varname, variable_list, freq='1Y', suffix='year_max')

class Reporting(object):
    def __init__(
            self,
            model,
            variable_list,
            num_samples=1,
            config_section='REPORTING',
            reporting_options=None,  # not used currently
            run_id=None              # not used currently
    ):
        self.model = model
        self.variable_list = variable_list
        self.outdir = self.model.config.output_directory
        self.reporting_variables = _get_reporting_variables(
            self.model.config, config_section
        )
        self.summary_variables = _get_summary_variables(
            self.model.config, 'MONTE_CARLO'
        )        
        self.num_samples = num_samples
        self.create_reporting_variables()

    def create_reporting_variables(self):
        self.output_variables = {}
        for sample in range(1, self.num_samples + 1):
            output_variables_dict = {}
            for option, varnames in self.reporting_variables.items():
                if varnames is not None:
                    for varname in varnames:
                        output_variables_dict[varname + '_' + option] = open_reporting_variable(
                            self.model,
                            varname,
                            option,
                            sample,
                            self.variable_list
                        )
                        
            self.output_variables[sample] = output_variables_dict
            
    def initial(self, sample=1):
        for _, value in self.output_variables[sample].items():
            value.initial()

    def dynamic(self, sample=1):
        for _, value in self.output_variables[sample].items():
            value.update()
            
    def create_mc_summary_variable(self):
        print("Hello, world")
        print(self.summary_variables)
        for option, varnames in self.summary_variables.items():
            if varnames is not None:
                for varname in varnames:
                    da_list = []
                    for sample in range(1, self.num_samples + 1):
                        obj = self.output_variables[sample][str(varname) + '_' + str(option)]                        
                        da_list.append(xarray.open_dataset(obj.filename)[obj._netcdf.attr.shortname])
                    da = xarray.concat(da_list, dim='run')
                    da_mean = xarr_merged.mean(dim='run')
                    da_std = xarr_merged.std(dim='run')
                    da_var = xarr_merged.var(dim='run')
                    # da.to_netcdf()
                    
