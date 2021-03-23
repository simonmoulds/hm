#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import sys
import numpy as np
import pickle

# These two functions are copied from LISFLOOD:
def dumpObject(name, var, num):
    path1 = os.path.join(str(num), 'stateVar', name)
    try:
        with open(path1, 'wb') as f:
            pickle.dump(var, f)
    except:
        pass

def loadObject(name, num):
    path1 = os.path.join(str(num), 'stateVar', name)
    try:
        with open(path1, 'rb') as f:
            return pickle.load(f)
    except FileNotFoundError:
        return 0
    
class stateVar(object):

    def __init__(self, dynamicmodel):
        self.dynamicmodel = dynamicmodel
        
    def initial(self):
        try:
            self.state_varnames = self.dynamicmodel.model.state_varnames
        except:
            self.state_varnames = []
        
    def dynamic(self):
        if self.dynamicmodel.currentTimeStep() in self.dynamicmodel.dump_timesteps:
            sample = str(self.dynamicmodel.currentSampleNumber())
            for varname in self.state_varnames:
                dumpObject(varname, vars(self.dynamicmodel.model)[varname], sample)                

    def resume(self):
        sample = str(self.dynamicmodel.currentSampleNumber())
        for varname in self.state_varnames:
            vars(self.dynamicmodel.model)[varname] = loadObject(varname, sample)
                
    
