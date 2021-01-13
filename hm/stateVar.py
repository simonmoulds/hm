#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import sys
import numpy as np
import pickle

# These two functions are copied from LISFLOOD:
def dumpObject(name, var, num):
    path1 = os.path.join(str(num), 'stateVar', name)
    file_object1 = open(path1, 'w')
    pickle.dump(var, file_object1)
    file_object1.close()

def loadObject(name, num):
    path1 = os.path.join(str(num), 'stateVar', name)
    file_handler1 = open(path1, 'w')
    var = pickle.load(file_handler1)
    file_handler1.close()
    return(var)
    
class stateVar(object):

    def __init__(self, dynamicmodel, timesteps):
        # timesteps : timesteps to r/w dump
        self.dynamicmodel = dynamicmodel
        self.timesteps = timesteps
        self.state_vars = dynamicmodel.state_vars

    def dynamic(self):
        # filter_timesteps = self.dynamicmodel.config.KALMAN_FILTER['filter_timesteps']
        if self.model.apply_kalman_filter and self.model.currentTimeStep() in self.timesteps:
            sample = str(self.currentSampleNumber())
            for varname in self.state_vars:
                try:
                    dumpObject(varname, vars(self.model)[varname], sample)
                except:
                    pass

    def resume(self):
        sample = str(self.currentSampleNumber())
        for varname in self.state_vars:
            try:
                loadObject(varname, sample)
            except:
                pass
                
    
