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

    def __init__(self, dynamicmodel):
        self.dynamicmodel = dynamicmodel
        try:
            self.state_varnames = dynamicmodel.model.state_varnames
        except:
            self.state_varnames = []

    def dynamic(self):
        # if self.model.apply_kalman_filter and self.model.currentTimeStep() in self.dynamicmodel.dump_timesteps:
        if self.dynamicmodel.currentTimeStep() in self.dynamicmodel.dump_timesteps:
            sample = str(self.dynamicmodel.currentSampleNumber())
            for varname in self.state_varnames:
                print('varname:', varname)
                print('value  :', vars(self.dynamicmodel.model)[varname])
                # try:
                dumpObject(varname, vars(self.dynamicmodel.model)[varname], sample)
                # except:
                #     pass

    def resume(self):
        sample = str(self.dynamicmodel.currentSampleNumber())
        for varname in self.state_varnames:
            # try:
            vars(self.dynamicmodel.model)[varname] = loadObject(varname, sample)
            print(varname, vars(self.dynamicmodel.model)[varname])
            # except:
            #     pass
                
    
