#!/usr/bin/env python
# -*- coding: utf-8 -*-

import time
import datetime
import netCDF4

import logging
logger = logging.getLogger(__name__)

class ModelTime(object):

    def __init__(self):
        object.__init__(self)
      
    #FIXME: use __init__
    def getStartEndTimeSteps(self,strStartTime,strEndTime,showNumberOfTimeSteps=True):
        # get startTime, endTime, nrOfTimeSteps 
        sd = str(strStartTime).split('-')
        self._startTime = datetime.date(int(sd[0]), int(sd[1]), int(sd[2]))
        ed = str(strEndTime).split('-')
        self._endTime = datetime.date(int(ed[0]), int(ed[1]), int(ed[2]))
        self._nrOfTimeSteps = 1 + (self.endTime - self.startTime).days
        # self._spinUpStatus = False
        if showNumberOfTimeSteps == True: logger.info("number of time steps: "+str(self._nrOfTimeSteps))
        self._monthIdx = 0 # monthly indexes since the simulation starts
        self._annuaIdx = 0 #  yearly indexes since the simulation starts

    #FIXME: use __init__
    # def getStartEndTimeStepsForSpinUp(self,strStartTime,noSpinUp,maxSpinUps):
    #     # get startTime, endTime, nrOfTimeSteps for SpinUps
    #     sd = str(strStartTime).split('-')
    #     self._startTime = datetime.date(int(sd[0]), int(sd[1]), int(sd[2]))
    #     self._endTime = datetime.date(int(sd[0])+1, int(sd[1]), int(sd[2])) -\
    #                    datetime.timedelta(days=1)
    #     self._nrOfTimeSteps = 1 + (self.endTime - self.startTime).days
    #     self._spinUpStatus = True
    #     self._noSpinUp   = noSpinUp
    #     self._maxSpinUps = maxSpinUps
    #     self._monthIdx = 0 # monthly indexes since the simulation starts
    #     self._annuaIdx = 0 #  yearly indexes since the simulation starts

    def setStartTime(self, date):
        self._startTime = date
        self._nrOfTimeSteps = 1 + (self.endTime - self.startTime).days

    def setEndTime(self, date):
        self._endTime = date
        self._nrOfTimeSteps = 1 + (self.endTime - self.startTime).days

    # @property
    # def spinUpStatus(self):
    #     return self._spinUpStatus

    @property    
    def startTime(self):
        return self._startTime
    
    @property    
    def endTime(self):
        return self._endTime

    @property    
    def currTime(self):
        return self._currTime

    @property    
    def day(self):
        return self._currTime.day

    @property    
    def doy(self):
        return self._currTime.timetuple().tm_yday

    @property
    def startTimeDOY(self):
        return self._startTime.timetuple().tm_yday
    
    @property    
    def month(self):
        return self._currTime.month
    
    @property    
    def year(self):
        return self._currTime.year

    @property    
    def timeStepPCR(self):
        return self._timeStepPCR
    
    @property    
    def monthIdx(self):
        return self._monthIdx
    
    @property    
    def annuaIdx(self):
        return self._annuaIdx
    
    @property    
    def nrOfTimeSteps(self):
        return self._nrOfTimeSteps
    
    @property
    def fulldate(self):
        return self._fulldate
    
    @property
    def isLeapYear(self):
        return self.year % 4 == 0 and (self.year % 100 != 0 or self.year % 400 == 0)

    @property
    def endTimeNum(self):
        return netCDF4.date2num(datetime.datetime(self._endTime.year, self._endTime.month, self._endTime.day), units="days since 1900-01-01 00:00:00")

    @property
    def currentYearStartNum(self):
        return netCDF4.date2num(datetime.datetime(self.year, 1, 1), units="days since 1900-01-01 00:00:00")

    def update(self,timeStepPCR):
        self._timeStepPCR = timeStepPCR
        self._currTime = self._startTime + datetime.timedelta(days=1 * (timeStepPCR - 1))
        
        #~ self._fulldate = str(self.currTime.strftime('%Y-%m-%d'))     # This does not work for the date before 1900
        self._fulldate = '%04i-%02i-%02i' %(self._currTime.year, self._currTime.month, self._currTime.day)
        
        # if self.spinUpStatus == True : 
        #     logger.info("Spin-Up "+str(self._noSpinUp)+" of "+str(self._maxSpinUps))

        # The following contains hours, minutes, seconds, etc. 
        self._currTimeFull = datetime.datetime(self.year,self.month,self.day)
        
        # check if a certain day is the last day of the month
        if self.isLastDayOfMonth():
            self._monthIdx = self._monthIdx + 1

        # check if a certain day is the last day of the year
        if self.isLastDayOfYear():
            self._annuaIdx = self._annuaIdx + 1
            
    def isFirstTimestep(self):
        return self.timeStepPCR == 1

    def isFirstDayOfMonth(self):
        return self.day == 1
    
    def isFirstDayOfYear(self):
        return self.doy== 1
    
    def isLastDayOfMonth(self):
        tomorrow = self.currTime + datetime.timedelta(days=1)
        
        #tomorrow is the first day of the month
        return tomorrow.day == 1
    
    def isLastDayOfYear(self):
        tomorrow = self.currTime + datetime.timedelta(days=1)
        
        #tomorrow is the first day of the year
        return tomorrow.timetuple().tm_yday == 1

    def isLastTimeStep(self):
        return self._currTime == self._endTime

    def yesterday(self):
        yesterday = self.currTime - datetime.timedelta(days=1)
        return str(yesterday.strftime('%Y-%m-%d'))

    #FIXME: use isLastDayOfMonth
    @property
    def endMonth(self):
        return self.isLastDayOfMonth()

    #FIXME: use isLastDayOfYear
    @property
    def endYear(self):
        return self.isLastDayOfYear()
    
    def __str__(self):
        return str(self._currTime)

def myfun(x: int) -> int: 
    return x * 10

class A(object):
    def __init__(self, x: int):
        self.x = x

class ModelTimeNew(object):

    def __init__(self, start_time, end_time, timedelta):
        # if isinstance(start_time, datetime.datetime):
        self._start_time = start_time
        self._end_time = end_time
        self._dt = timedelta
        if (end_time - start_time) % timedelta == datetime.timedelta(0):
            self._n_timesteps = (end_time - start_time) / timedelta
        else:
            msg = "Given end_time is not a multiple of timedelta"
            raise ValueError(msg)
        
        self._times = pd.date_range(start_time, periods=self._n_timesteps).to_pydatetime().tolist()
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
    def start_time(self):
        return self._start_time
    
    @property    
    def end_time(self):
        return self._end_time

    @property
    def dt(self):
        return self._dt
    
    @property    
    def curr_time(self):
        return self._curr_time

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
        # self._curr_time = self._start_time + datetime.timedelta(days=1 * (timeStepPCR - 1))
        self._curr_time = self._times[self._timestep]        
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
        return self._curr_time == self._end_time

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
        return self._n_timesteps
