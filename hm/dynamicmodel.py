#!/usr/bin/env python
# -*- coding: utf-8 -*-

from .reporting import Reporting, DummyReporting
from .pcraster.dynamicPCRasterBase import DynamicModel
from .pcraster.mcPCRasterBase import MonteCarloModel
from .pcraster.kfPCRasterBase import EnKfModel

class HmDynamicModel(DynamicModel):    
    def __init__(
        self,
        model,
        config,
        modeltime,
        domain,
        variable_list,
        init = None
    ):
        DynamicModel.__init__(self)
        self.modeltime = modeltime
        self.model = model(
            config,
            modeltime,
            domain,
            init
        )
        self.model.initial()
        if config.REPORTING['report'] == True:
            self.reporting = Reporting(self.model, variable_list)
        else:
            self.reporting = DummyReporting()

    def initial(self):
        self.reporting.initial()

    def dynamic(self):
        self.model.time.update(self.currentTimeStep())
        self.model.dynamic()
        self.reporting.dynamic()

class HmMonteCarloModel(DynamicModel, MonteCarloModel):
    def __init__(
            self,
            model,
            config,
            modeltime,
            domain,
            variable_list,
            init = None
    ):
        DynamicModel.__init__(self)
        MonteCarloModel.__init__(self)
        self.modeltime = modeltime
        self.model = model(
            config,
            modeltime,
            domain,
            init
        )
        self.model.initial()        
        self.config = config
        self.variable_list = variable_list

    def premcloop(self):        
        if self.config.REPORTING['report'] == True:
            self.reporting = Reporting(self.model, self.variable_list, self.nrSamples())
        else:
            self.reporting = DummyReporting()
        
    def initial(self):        
        self.reporting.initial(self.currentSampleNumber())

    def dynamic(self):
        self.model.time.update(self.currentTimeStep())
        self.model.dynamic()
        self.reporting.dynamic(self.currentSampleNumber())

    def postmcloop(self):
        pass
        
class HmEnKfModel(DynamicModel, MonteCarloModel, EnKfModel):
    def __init__(
            self,
            model,
            config,
            modeltime,
            domain,
            variable_list,
            init = None
    ):
        DynamicModel.__init__(self)
        MonteCarloModel.__init__(self)
        EnKfModel.__init__(self)

    def premcloop(self):
        pass
    
    def initial(self):
        pass

    def dynamic(self):
        pass

    def postmcloop(self):
        pass

# from .reporting import Reporting, DummyReporting
# import logging
# logger = logging.getLogger(__name__)

# class DynamicModel(object):
#   def __init__(self):
#     self.silentModelOutput = False
#     self._d_nrTimeSteps = 0
#     self.currentStep = 0
#     self._d_firstTimeStep = 1
#     self.inTimeStep = False
#     self.inInitial = False
#     self.inDynamic = False
	
#   def initial(self):
#     #print "Implement 'initial' method"
#     ii =1
	
#   def dynamic(self):
#     print("Implement 'dynamic' method")

#   def timeSteps(self):
#     """
#     Return a list of time steps
#     """
#     return list(range(self.firstTimeStep(), self.nrTimeSteps() + 1))

#   def nrTimeSteps(self):
#     """
#     Return the number of time steps
#     """
#     assert self._d_nrTimeSteps
#     return self._d_nrTimeSteps

#   def currentTimeStep(self):
#     """
#     Return the current time step in the range from firstTimeStep to nrTimeSteps.
#     """
#     assert self.currentStep >= 0
#     return self.currentStep

#   def firstTimeStep(self):
#     """
#     Return first timestep of a model.
#     """
#     assert self._d_firstTimeStep
#     return self._d_firstTimeStep

#   def setQuiet(self,
#     quiet=True):
#     """
#     Disables the progress display of timesteps.
#     """
#     self.silentModelOutput = quiet	
	
#   def _silentModelOutput(self):
#     return self.silentModelOutput

#   def _inDynamic(self):
#     return self.inDynamic

#   def _inInitial(self):
#     return self.inInitial

#   def _setInInitial(self, value):
#     assert isinstance(value, bool)
#     self.inInitial = value

#   def _setInDynamic(self, value):
#     assert isinstance(value, bool)
#     self.inDynamic = value

#   def _inTimeStep(self):
#     """
#     Returns whether a time step is currently executing.
#     """
#     #if hasattr(self._userModel(), "_d_inTimeStep"):
#     #  return self._userModel()._d_inTimeStep
#     #else:
#     #  return False
#     return self.inTimeStep

#   def _setInTimeStep(self, value):
#     assert isinstance(value, bool)
#     self.inTimeStep = value	
  
#   def _setFirstTimeStep(self, firstTimeStep):

#     if not isinstance(firstTimeStep, int):
#       msg = "first timestep argument (%s) of DynamicFramework must be of type int" % (type(firstTimeStep))
#       raise AttributeError(msg)

#     if firstTimeStep <= 0:
#       msg = "first timestep argument (%s) of DynamicFramework must be > 0" % (firstTimeStep)
#       raise AttributeError(msg)

#     if firstTimeStep > self.nrTimeSteps():
#       msg = "first timestep argument (%s) of DynamicFramework must be smaller than given last timestep (%s)" % (firstTimeStep, self.nrTimeSteps())
#       raise AttributeError(msg)

#     self._d_firstTimeStep = firstTimeStep
	
#   def _setNrTimeSteps(self, TimeStep):
#     """
#     Set the number of time steps to run.
#     """
#     if not isinstance(TimeStep, int):
#       msg = "last timestep argument (%s) of DynamicFramework must be of type int" % (type(TimeStep))
#       raise AttributeError(msg)

#     if TimeStep <= 0:
#       msg = "last timestep argument (%s) of DynamicFramework must be > 0" % (TimeStep)
#       raise AttributeError(msg)

#     self._d_nrTimeSteps = TimeStep

#   def _setCurrentTimeStep(self,
#     step):
#     """
#     Set the current time step.
#     """
#     if step <= 0:
#       msg = "Current timestep must be > 0"
#       raise AttributeError(msg)

#     if step > self.nrTimeSteps():
#       msg = "Current timestep must be <= %d (nrTimeSteps)"
#       raise AttributeError(msg)

#     self.currentStep = step

# class HmDynamicModel(DynamicModel):    
#     def __init__(
#         self,
#         model,
#         config,
#         modeltime,
#         domain,
#         variable_list,
#         init = None
#     ):
#         """Class for running a hydrological model.
        
#         Parameters
#         ----------
#         model : class
#         config : 
#         modeltime : ModelTime 
#         init : 
#         report : boolean
#         """
#         DynamicModel.__init__(self)
#         self.modeltime = modeltime
#         self.model = model(
#             config,
#             modeltime,
#             domain,
#             init
#         )
#         self.model.initial()
#         if config.REPORTING['report'] == 'True':
#             self.reporting = Reporting(self.model, variable_list)
#         else:
#             self.reporting = DummyReporting()

#     def initial(self):
#         self.reporting.initial()

#     def dynamic(self):
#         self.model.time.update(self.currentTimeStep())
#         self.model.dynamic()
#         # self.reporting.report()
#         self.reporting.dynamic()
