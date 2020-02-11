#!/usr/bin/env python
# -*- coding: utf-8 -*-

import time
import datetime
import pandas as pd

import logging
logger = logging.getLogger(__name__)

class ModelTime(object):

    def __init__(
            self,
            starttime,
            endtime,
            timedelta
    ):
        """To create... 

        Parameters
        ----------
        starttime : pandas.Timestamp
            The start time of the model simulation
        endtime : pandas.Timestamp
            The end time of the simulation
        timedelta : pandas.Timedelta
            The model time step
        """
        self._starttime = pd.Timestamp(starttime)
        self._endtime = pd.Timestamp(endtime)
        self._dt = pd.Timedelta(timedelta)
        if (endtime - starttime) % timedelta == datetime.timedelta(0):
            self._n_timesteps = (endtime - starttime) / timedelta
        else:
            raise ValueError(
                'Given endtime is not a multiple of timedelta'
            )
        self._times = pd.date_range(starttime, endtime, periods=self._n_timesteps + 1)
        # if show_number_of_timesteps:
        #     logger.info("number of time steps: " + str(self._n_timesteps))
        self._timestep = 0
        self._month_index = 0
        self._year_index = 0

    # def setStartTime(self, date):
    #     self._startTime = date
    #     self._nrOfTimeSteps = 1 + (self.endTime - self.startTime).days

    # def setEndTime(self, date):
    #     self._endTime = date
    #     self._nrOfTimeSteps = 1 + (self.endTime - self.startTime).days

    # @property
    # def spinUpStatus(self):
    #     return self._spinUpStatus

    @property    
    def starttime(self):
        return self._starttime
    
    @property    
    def endtime(self):
        return self._endtime

    @property
    def dt(self):
        return self._dt

    @property
    def times(self):
        return self._times
    @property    
    def curr_time(self):
        return self._curr_time
    @property
    def timestamp(self):
        return pd.Timestamp(self._curr_time)
    @property    
    def day(self):
        return self._curr_time.day

    @property    
    def doy(self):
        return self._curr_time.timetuple().tm_yday

    # @property
    # def startTimeDOY(self):
    #     return self._startTime.timetuple().tm_yday
    
    @property    
    def month(self):
        return self._curr_time.month
    
    @property    
    def year(self):
        return self._curr_time.year

    # @property    
    # def timeStepPCR(self):
    #     return self._timeStepPCR
    @property
    def timestep(self):
        return self._timestep
    
    # @property    
    # def monthIdx(self):
    #     return self._monthIdx
    
    # @property    
    # def annuaIdx(self):
    #     return self._annuaIdx
    
    # @property    
    # def nrOfTimeSteps(self):
    #     return self._nrOfTimeSteps
    
    @property
    def fulldate(self):
        return self._fulldate
    
    @property
    def is_leap_year(self):
        return self.year % 4 == 0 and (self.year % 100 != 0 or self.year % 400 == 0)

    # @property
    # def endTimeNum(self):
    #     return netCDF4.date2num(datetime.datetime(self._endTime.year, self._endTime.month, self._endTime.day), units="days since 1900-01-01 00:00:00")

    # @property
    # def currentYearStartNum(self):
    #     return netCDF4.date2num(datetime.datetime(self.year, 1, 1), units="days since 1900-01-01 00:00:00")

    def update(self, timestep):
        self._timestep = timestep
        # self._curr_time = self._starttime + datetime.timedelta(days=1 * (timeStepPCR - 1))
        self._curr_time = self._times[(self._timestep - 1)]        
        # self._fulldate = '%04i-%02i-%02i' % (self._curr_time.year, self._curr_time.month, self._curr_time.day)
        self._fulldate = self._curr_time.strftime("%Y-%m-%d %H:%M:%S")
        # if self.spinUpStatus == True : 
        #     logger.info("Spin-Up "+str(self._noSpinUp)+" of "+str(self._maxSpinUps))
        # The following contains hours, minutes, seconds, etc. 
        # self._currTimeFull = datetime.datetime(self.year,self.month,self.day)

        # TODO: this should be last timestep of month/year
        if self.is_last_day_of_month():
            self._month_index = self._month_index + 1
        if self.is_last_day_of_year():
            self._year_index = self._year_index + 1

    @property
    def is_first_timestep(self):
        return self.timestep == 1

    def is_first_day_of_month(self):
        return self.day == 1
    
    def is_first_day_of_year(self):
        return self.doy == 1
    
    def is_last_day_of_month(self):
        tomorrow = self.curr_time + datetime.timedelta(days=1)
        return tomorrow.day == 1
    
    def is_last_day_of_year(self):
        tomorrow = self.curr_time + datetime.timedelta(days=1)
        return tomorrow.timetuple().tm_yday == 1

    def is_last_timestep(self):
        return self._curr_time == self._endtime

    def yesterday(self):
        yesterday = self.curr_time - datetime.timedelta(days=1)
        return str(yesterday.strftime('%Y-%m-%d'))

    # #FIXME: use isLastDayOfMonth
    # @property
    # def endMonth(self):
    #     return self.isLastDayOfMonth()

    # #FIXME: use isLastDayOfYear
    # @property
    # def endYear(self):
    #     return self.isLastDayOfYear()
    
    def __str__(self):
        return str(self._curr_time)

    def __len__(self):
        return int(self._n_timesteps)
