#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import shutil
import numpy as np
import datetime

from box import Box
# from collections import OrderedDict
# from .OutputNetCDF import OutputNetCDF
# from . import file_handling 
from .constants import *
from .utils import *

import logging
logger = logging.getLogger(__name__)

def get_variable_names_for_reporting(config, option):
    try:
        return [str(var.strip()) for var in config.REPORTING[option].split(',')]
    except KeyError:
        return []

def get_reporting_variables(config):
    # var_dict = {'all' : []}
    var_dict = {}
    all_vars = []
    for option in allowed_reporting_options:
        var_names = get_variable_names_for_reporting(config, option)
        var_dict[option] = var_names
        all_vars += var_names
    all_vars = sorted(set(all_vars))
    var_dict['all'] = all_vars
    return Box(var_dict, frozen_box=True)

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
    return Box(attr, frozen_box=True)

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
        self.data += vars(self._model)[self.varname]
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
    
class reporting(object):
    def __init__(
            self,
            model,
            output_dir='.',
            netcdf_attr=None
            # reporting_options=None,
            # run_id=None
    ):
        self._model = model
        self._outdir = output_dir
        self._reporting_variables = get_reporting_variables(self._model.config)
        self.initiate_reporting()

    def create_output_variables(self):
        """Create arrays for summary outputs."""
        for option, varnames in self._reporting_variables.items():
            if varnames is not None:
                for var in varnames:
                    # create array for summary output
                    var_shape = self._model[var].shape[:-1] + (self._model.domain.mask.shape)
                    vars(self)[var + '_' + option] = np.ma.array(
                        np.zeros(var_shape),
                        mask=np.broadcast_to(self._model.domain.mask, arr.shape)
                    )
                    self.create_netcdf(var)

    def create_netcdf(self):
        # nc.Dataset(filename, 'w', format)
        pass
        
    def report(self):
        logger.info("reporting for time %s", self._modelTime.currTime)
        self.post_processing()
        # self.time_stamp = datetime.datetime(
        #     self._modelTime.year,
        #     self._modelTime.month,
        #     self._modelTime.day,
        #     0)        
        self.report_daily_total()
        self.report_month_average()
        self.report_month_end()
        self.report_month_total()
        self.report_month_maximum()
        self.report_year_average()
        self.report_year_end()
        self.report_year_total()
        self.report_year_maximum()
        
    # def post_processing(self):
    #     """Copy model variables."""
    #     # TODO: I think this class is redundant: REMOVE
    #     if len(self.variables_for_report) > 0:
    #         for var in self.variables_for_report.all:
    #             d = vars(self._model)[var].copy()
    #             # create an array with correct spatial dimension
    #             arr = np.zeros(
    #                 d.shape[:-1] + (self._model.domain.mask.shape)
    #             )
    #             arr[..., self._model.domain.mask] = d
    #             # promote array to masked array; this will detect
    #             # the default fill value based on arr.dtype
    #             arr = np.ma.array(
    #                 arr,
    #                 mask=np.broadcast_to(self._model.domain.mask, arr.shape)
    #             )
    #             vars(self)[var] = arr.copy()

    # TODO: make these part of modeltime object
    # def get_n_timestep_in_period(self):
    #     n_timestep_in_three_hourly_period
    #     n_timestep_in_hourly_period
    #     n_timestep_in_daily_period
    #     n_timestep_in_weekly_period
    #     n_timestep_in_dekadal_period
    #     n_timestep_in_monthly_period
    #     n_timestep_in_annual_period
        
    def report(self, period, is_end_time, n_timesteps_in_period):
        n_timesteps_in_period = np.min(n_timesteps_in_period, self._model.model_time.timeStepPCR)
        for varname in self._reporting_variables.all:
            var = reshape_array(vars(self._model)[varname], self._model.domain.model)
            if varname in self._
        for varname in self._reporting_variables[period + '_mean']:
            self.report_mean(
                var,
                varname + '_' + period + '_mean',
                is_end_time,
                n_timesteps_in_period
            )            
        for varname in self._reporting_variables[period + '_min']:
            var = reshape_array(vars(self._model)[varname], self._model.domain.model)
            self.report_min(
                varname,
                varname + '_' + period + '_min',
                is_end_time,
                n_timesteps_in_period
            )                
        for varname in self._reporting_variables[period + '_max']:
            var = reshape_array(vars(self._model)[varname], self._model.domain.model)
            self.report_max(
                varname,
                varname + '_' + period + '_max',
                is_end_time,
                n_timesteps_in_period
            )
        for varname in self._reporting_variables[period + '_end']:
            var = reshape_array(vars(self._model)[varname], self._model.domain.model)
            self.report_end(
                varname,
                varname + '_' + period + '_end',
                is_end_time,
                n_timesteps_in_period
            )
        for varname in self._reporting_variables[period + '_total']:
            var = reshape_array(vars(self._model)[varname], self._model.domain.model)
            self.report_end(
                varname,
                varname + '_' + period + '_total',
                is_end_time,
                n_timesteps_in_period
            )
    
    def report_mean(self, var, varname_mean, is_end_period, n_timestep):
        vars(self)[varname_mean] += var
        if is_end_period:
            vars(self)[varname_mean] = vars(self)[varname_mean] / n_timestep
            # filename
            # add data to netcdf
            # reset
            
    def report_min(self, var, varname_min, is_end_period, n_timestep):
        vars(self)[varname_min] = vars(self)[varname_min].clip(None, var)
        if is_end_period:
            # filename
            # add data to netcdf
            # reset
            pass
    def report_max(self, var, varname_max, is_end_period, n_timestep):
        vars(self)[varname_max] = vars(self)[varname_max].clip(var, None)
        if is_end_period:
            # filename
            # add data to netcdf
            # reset
            pass
    def report_end(self, var, varname_end, is_end_period, n_timestep):
        if is_end_period:
            vars(self)[varname_end] = var
            # filename
            # add data to netcdf
            # reset
    def report_total(self, var, varname_max, is_end_period, n_timestep):
        vars(self)[varname_total] += var
        if is_end_period:
            # filename
            # add data to netcdf
            # reset
            pass
            
# class Reporting(object):

#     def __init__(self, model, output_dir, netcdf_attr, reporting_options, variable_list, run_id=None):
#         self._model = model
#         self._modelTime = model._modelTime
#         self.output_dir = output_dir
#         self.reporting_options = reporting_options
#         self.initiate_reporting(netcdf_attr, variable_list, run_id)        

#     def create_netcdf_file(self, var, suffix):
#         ncFile = self.output_dir + "/" + str(var) + str(suffix) + ".nc"
#         self.netcdfObj.create_netCDF(ncFile, var)

#     def initiate_reporting(self, netcdf_attr, variable_list, run_id):
#         """Function to create netCDF files for each output 
#         variable
#         """
#         self.netcdfObj = OutputNetCDF(
#             netcdf_attr,
#             self._model.dimensions,
#             variable_list)

#         if run_id is None:
#             run_id = ''
#         else:
#             run_id = '_' + str(run_id)
#         self.run_id = run_id

#         self.initiate_daily_total_reporting()
#         self.initiate_month_average_reporting()
#         self.initiate_month_end_reporting()
#         self.initiate_month_total_reporting()
#         self.initiate_month_maximum_reporting()
#         self.initiate_year_average_reporting()
#         self.initiate_year_end_reporting()
#         self.initiate_year_total_reporting()
#         self.initiate_year_maximum_reporting()
        
#         self.variables_for_report = (
#             self.outDailyTotal +
#             self.outMonthAvgNC +
#             self.outMonthEndNC +
#             self.outMonthTotNC +
#             self.outMonthMaxNC +
#             self.outYearAvgNC +
#             self.outYearEndNC +
#             self.outYearTotNC +
#             self.outYearMaxNC
#             )

#         # reduce above list to unique values, and remove None
#         self.variables_for_report = list(set(self.variables_for_report))
#         if "None" in self.variables_for_report:
#             self.variables_for_report.remove("None")

#     def get_variable_names_for_reporting(self, option):
#         var_names = [str(var.strip()) for var in self.reporting_options[option].split(',')]
#         return list(set(var_names))
    
#     def initiate_daily_total_reporting(self):
#         self.outDailyTotal = ["None"]
#         try:
#             self.outDailyTotal = self.get_variable_names_for_reporting('outDailyTotal')
#         except:
#             pass
        
#         if self.outDailyTotal[0] != "None":
#             for var in self.outDailyTotal:
#                 logger.info("Creating the netcdf file for reporting the daily value of variable %s.", str(var))
#                 self.create_netcdf_file(var, self.run_id + "_dailyTot_output")

#     def initiate_month_average_reporting(self):
#         self.outMonthAvgNC = ["None"]
#         try:
#             self.outMonthAvgNC = self.get_variable_names_for_reporting('outMonthAvgNC')
#         except:
#             pass
#         if self.outMonthAvgNC[0] != "None":
#             for var in self.outMonthAvgNC:
#                 logger.info("Creating the netcdf file for reporting the monthly average of variable %s.", str(var))
#                 self.create_netcdf_file(var, self.run_id + "_monthAvg_output")
#                 vars(self)[var+'_monthAvg'] = np.zeros(vars(self._model)[var].shape[:-1] + (self._model.nLat, self._model.nLon))

#     def initiate_month_end_reporting(self):
#         self.outMonthEndNC = ["None"]
#         try:
#             self.outMonthEndNC = list(set(self.reporting_options['outMonthEndNC'].split(",")))
#         except:
#             pass
#         if self.outMonthEndNC[0] != "None":
#             for var in self.outMonthEndNC:
#                 logger.info("Creating the netcdf file for reporting the month end value of variable %s.", str(var))
#                 self.create_netcdf_file(var, self.run_id + "_monthEnd_output")
#                 vars(self)[var+'_monthEnd'] = np.zeros(vars(self._model)[var].shape[:-1] + (self._model.nLat, self._model.nLon))

#     def initiate_month_total_reporting(self):
#         self.outMonthTotNC = ["None"]
#         try:
#             self.outMonthTotNC = self.get_variable_names_for_reporting('outMonthTotNC')
#         except:
#             pass
#         if self.outMonthTotNC[0] != "None":
#             for var in self.outMonthTotNC:
#                 logger.info("Creating the netcdf file for reporting the monthly average of variable %s.", str(var))
#                 self.create_netcdf_file(var, self.run_id + "_monthTot_output")
#                 vars(self)[var+'_monthTot'] = np.zeros(vars(self._model)[var].shape[:-1] + (self._model.nLat, self._model.nLon))

#     def initiate_month_maximum_reporting(self):
#         self.outMonthMaxNC = ["None"]
#         try:
#             self.outMonthMaxNC = list(set(self.reporting_options['outMonthMaxNC'].split(",")))
#         except:
#             pass
#         if self.outMonthMaxNC[0] != "None":
#             for var in self.outMonthMaxNC:
#                 logger.info("Creating the netcdf file for reporting the monthly maximum of variable %s.", str(var))
#                 self.create_netcdf_file(var, self.run_id + "_monthMax_output")
#                 vars(self)[var+'_monthMax'] = np.zeros(vars(self._model)[var].shape[:-1] + (self._model.nLat, self._model.nLon))

#     def initiate_year_average_reporting(self):
#         self.outYearAvgNC = ["None"]
#         try:
#             self.outYearAvgNC = list(set(self.reporting_options['outYearAvgNC'].split(",")))
#         except:
#             pass
#         if self.outYearAvgNC[0] != "None":
#             for var in self.outYearAvgNC:
#                 logger.info("Creating the netcdf file for reporting the yearly average of variable %s.", str(var))
#                 self.create_netcdf_file(var, self.run_id + "_yearAvg_output")
#                 vars(self)[var+'_yearAvg'] = np.zeros(vars(self._model)[var].shape[:-1] + (self._model.nLat, self._model.nLon))
#                 # vars(self)[var+'_yearAvg'] = np.zeros((vars(self._model)[var].shape))

#     def initiate_year_end_reporting(self):
#         self.outYearEndNC = ["None"]
#         try:
#             self.outYearEndNC = list(set(self.reporting_options['outYearEndNC'].split(",")))
#         except:
#             pass
#         if self.outYearEndNC[0] != "None":
#             for var in self.outYearEndNC:
#                 logger.info("Creating the netcdf file for reporting the year end value of variable %s.", str(var))
#                 self.create_netcdf_file(var, self.run_id + "_yearEnd_output")
#                 vars(self)[var+'_yearEnd'] = np.zeros((vars(self._model)[var].shape))

#     def initiate_year_total_reporting(self):
#         self.outYearTotNC = ["None"]
#         try:
#             self.outYearTotNC = self.get_variable_names_for_reporting('outYearTotNC')
#         except:
#             pass
#         if self.outYearTotNC[0] != "None":
#             for var in self.outYearTotNC:
#                 logger.info("Creating the netcdf file for reporting the yearly total of variable %s.", str(var))
#                 self.create_netcdf_file(var, self.run_id + "_yearTot_output")
#                 vars(self)[var+'_yearTot'] = np.zeros(vars(self._model)[var].shape[:-1] + (self._model.nLat, self._model.nLon))

#     def initiate_year_maximum_reporting(self):
#         self.outYearMaxNC = ["None"]
#         try:
#             self.outYearMaxNC = list(set(self.reporting_options['outYearMaxNC'].split(",")))
#         except:
#             pass
#         if self.outYearMaxNC[0] != "None":
#             for var in self.outYearMaxNC:
#                 logger.info("Creating the netcdf file for reporting the yearly maximum of variable %s.", str(var))
#                 self.create_netcdf_file(var, self.run_id + "_yearMax_output")
#                 vars(self)[var+'_yearMax'] = np.zeros(vars(self._model)[var].shape[:-1] + (self._model.nLat, self._model.nLon))
                
#     def post_processing(self):
#         """Function to process model variables to output variables. 
#         This generally means assigning values to an array with 
#         spatial dimensions equal to those of the current landmask.
#         """
#         if len(self.variables_for_report) > 0:
#             for var in self.variables_for_report:
#                 d = vars(self._model)[var].copy()
#                 arr = np.full(
#                     d.shape[:-1] + (self._model.nLat, self._model.nLon),
#                     file_handling.missing_value
#                     )
#                 arr[...,self._model.landmask] = d
#                 vars(self)[var] = arr.copy()

#     def report(self):
#         logger.info("reporting for time %s", self._modelTime.currTime)
#         self.post_processing()
#         self.time_stamp = datetime.datetime(
#             self._modelTime.year,
#             self._modelTime.month,
#             self._modelTime.day,
#             0)
        
#         self.report_daily_total()
#         self.report_month_average()
#         self.report_month_end()
#         self.report_month_total()
#         self.report_month_maximum()
#         self.report_year_average()
#         self.report_year_end()
#         self.report_year_total()
#         self.report_year_maximum()
        
#     def report_daily_total(self):
#         if self.outDailyTotal[0] != "None":
#             for var in self.outDailyTotal:
#                 fn = self.output_dir + "/" + str(var) + self.run_id + "_dailyTot_output.nc"
#                 self.netcdfObj.add_data_to_netcdf(
#                     fn,
#                     var,
#                     self.__getattribute__(var),
#                     self.time_stamp)

#     def report_month_average(self):
#         if self.outMonthAvgNC[0] != "None":
#             for var in self.outMonthAvgNC:
#                 if self._modelTime.day == 1: vars(self)[var+'_monthAvg'].fill(0)
#                 vars(self)[var+'_monthAvg'] += vars(self)[var]
#                 if self._modelTime.endMonth:
#                     divd = np.min((self._modelTime.timeStepPCR, self._modelTime.day))
#                     vars(self)[var+'_monthAvg'] = vars(self)[var+'_monthAvg'] / divd
#                     fn = self.output_dir + "/" + str(var) + self.run_id + "_monthAvg_output.nc"
#                     self.netcdfObj.add_data_to_netcdf(
#                         fn,
#                         var,
#                         self.__getattribute__(var+'_monthAvg'),
#                         self.time_stamp)
                    
#     def report_month_end(self):
#         if self.outMonthEndNC[0] != "None":
#             for var in self.outMonthEndNC:
#                 if self._modelTime.endMonth:
#                     vars(self)[var+'_monthEnd'] = vars(self)[var]
#                     fn = self.output_dir + "/" + str(var) + self.run_id + "_monthEnd_output.nc"
#                     self.netcdfObj.add_data_to_netcdf(
#                         fn,
#                         var,
#                         self.__getattribute__(var+'_monthEnd'),
#                         self.time_stamp)

#     def report_month_total(self):
#         if self.outMonthTotNC[0] != "None":
#             for var in self.outMonthTotNC:
#                 if self._modelTime.day == 1: vars(self)[var+'_monthTot'].fill(0)
#                 vars(self)[var+'_monthTot'][...,self._model.landmask] += vars(self)[var][...,self._model.landmask]
                
#                 if self._modelTime.endMonth:
#                     fn = self.output_dir + "/" + str(var) + self.run_id + "_monthTot_output.nc"
#                     self.netcdfObj.add_data_to_netcdf(
#                         fn,
#                         var,
#                         self.__getattribute__(var+'_monthTot'),
#                         self.time_stamp)                    

#     def report_month_maximum(self):
#         if self.outMonthMaxNC[0] != "None":
#             for var in self.outMonthMaxNC:
#                 if self._modelTime.day == 1: vars(self)[var+'_monthMax'].fill(0)
#                 vars(self)[var+'_monthMax'] = vars(self)[var+'_monthMax'].clip(vars(self)[var], None)
#                 if self._modelTime.endMonth:
#                     fn = self.output_dir + "/" + str(var) + self.run_id + "_monthMax_output.nc"
#                     self.netcdfObj.add_data_to_netcdf(
#                         fn,
#                         var,
#                         self.__getattribute__(var+'_monthMax'),
#                         self.time_stamp)                    

#     def report_year_average(self):
#         if self.outYearAvgNC[0] != "None":
#             for var in self.outYearAvgNC:
#                 if self._modelTime.doy == 1: vars(self)[var+'_yearAvg'].fill(0)
#                 vars(self)[var+'_yearAvg'] += vars(self)[var]
#                 if self._modelTime.endYear:
#                     divd = np.min((self._modelTime.timeStepPCR, self._modelTime.doy))
#                     vars(self)[var+'_yearAvg'] = vars(self)[var+'_yearAvg'] / divd
#                     fn = self.output_dir + "/" + str(var) + self.run_id + "_yearAvg_output.nc"
#                     self.netcdfObj.add_data_to_netcdf(
#                         fn,
#                         var,
#                         self.__getattribute__(var+'_yearAvg'),
#                         self.time_stamp)
                    

#     def report_year_end(self):
#         if self.outYearEndNC[0] != "None":
#             for var in self.outYearEndNC:
#                 if self._modelTime.endYear:
#                     vars(self)[var+'_yearEnd'] = vars(self)[var]
#                     fn = self.output_dir + "/" + str(var) + self.run_id + "_yearEnd_output.nc"
#                     self.netcdfObj.add_data_to_netcdf(

#                         fn,
#                         var,
#                         self.__getattribute__(var+'_yearEnd'),
#                         self.time_stamp)

#     def report_year_total(self):
#         if self.outYearTotNC[0] != "None":
#             for var in self.outYearTotNC:
#                 if self._modelTime.doy == 1: vars(self)[var+'_yearTot'].fill(0)
#                 vars(self)[var+'_yearTot'][...,self._model.landmask] += vars(self)[var][...,self._model.landmask]
#                 if self._modelTime.endYear:
#                     fn = self.output_dir + "/" + str(var) + self.run_id + "_yearTot_output.nc"
#                     self.netcdfObj.add_data_to_netcdf(
#                         fn,
#                         var,
#                         self.__getattribute__(var+'_yearTot'),
#                         self.time_stamp)                    

#     def report_year_maximum(self):
#         if self.outYearMaxNC[0] != "None":
#             for var in self.outYearMaxNC:
#                 if self._modelTime.doy == 1:
#                     vars(self)[var+'_yearMax'].fill(0)
#                 vars(self)[var+'_yearMax'] = vars(self)[var+'_yearMax'].clip(vars(self)[var], None)
#                 if self._modelTime.endYear or self._modelTime.isLastTimeStep():
#                     fn = self.output_dir + "/" + str(var) + self.run_id + "_yearMax_output.nc"
#                     self.netcdfObj.add_data_to_netcdf(
#                         fn,
#                         var,
#                         self.__getattribute__(var+'_yearMax'),
#                         self.time_stamp)
