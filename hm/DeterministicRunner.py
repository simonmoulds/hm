#!/usr/bin/env python
# -*- coding: utf-8 -*-

from .dynamicPCRasterBase import DynamicModel
from .Reporting import Reporting

import logging
logger = logging.getLogger(__name__)

class DummyReporting(object):
    def __init__(self):
        pass

    def report(self):
        pass
    
class DeterministicRunner(DynamicModel):    
    def __init__(self, model, configuration, modelTime, variable_list, initialState = None, suffix = None, report = True):
        DynamicModel.__init__(self)
        self.modelTime = modelTime
        self.model = model(
            configuration,
            modelTime,
            initialState)
        self.model.initial()

        if report:
            self.reporting = Reporting(
                self.model,
                configuration.outNCDir,
                configuration.NETCDF_ATTRIBUTES,
                configuration.REPORTING,
                variable_list,
                suffix)
        else:
            self.reporting = DummyReporting()

    def initial(self):
        pass

    def dynamic(self):
        self.modelTime.update(self.currentTimeStep())
        self.model.dynamic()
        self.reporting.report()
