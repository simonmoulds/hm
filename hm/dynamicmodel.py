#!/usr/bin/env python
# -*- coding: utf-8 -*-

from pcraster.dynamicPCRasterBase import DynamicModel
from .Reporting import Reporting, DummyReporting

import logging
logger = logging.getLogger(__name__)

class HmDynamicModel(DynamicModel):    
    def __init__(
            self,
            model,
            config,
            modeltime,
            init = None
    ):
        """Class for running a hydrological model.
        
        Parameters
        ----------
        model : class
        config : 
        modeltime : ModelTime 
        init : 
        report : boolean
        """
        DynamicModel.__init__(self)
        self.modeltime = modeltime
        self.model = model(
            config,
            modelTime,
            init)
        self.model.initial()
        if config.REPORTING['report']:
            self.reporting = Reporting(self.model)
        else:
            self.reporting = DummyReporting()

    def initial(self):
        pass

    def dynamic(self):
        self.modelTime.update(self.currentTimeStep())
        self.model.dynamic()
        self.reporting.report()
