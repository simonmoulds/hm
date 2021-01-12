#!/usr/bin/env python
# -*- coding: utf-8 -*-

# ######################################################## #
# constants for input netCDF files                         #
# ######################################################## #

missing_value = -9999
allowed_xy_dim_names = ['xy', 'space']
allowed_y_dim_names = ['lat', 'latitude', 'y']
allowed_x_dim_names = ['lon', 'longitude', 'x']
allowed_z_dim_names = ['z', 'depth']
allowed_t_dim_names = ['time', 'tstep']
controlled_dim_names = allowed_y_dim_names + \
    allowed_x_dim_names + allowed_z_dim_names + allowed_t_dim_names

# ######################################################## #
# constants for output netCDF files                        #
# ######################################################## #

intervals = ['hourly', 'three_hourly', 'daily', 'dekadal', 'month', 'year']
methods = ['mean', 'max', 'min', 'end', 'total']
allowed_reporting_options = []
allowed_summary_options = []
for interval in intervals:
    for method in methods:
        allowed_reporting_options.append(interval + '_' + method)
        allowed_summary_options.append(interval + '_' + method + '_summary')
