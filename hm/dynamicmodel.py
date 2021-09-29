#!/usr/bin/env python
# -*- coding: utf-8 -*-

from .reporting import Reporting, DummyReporting
from .stateVar import stateVar
# from .pcraster.dynamicPCRasterBase import DynamicModel
# from .pcraster.mcPCRasterBase import MonteCarloModel
# from .pcraster.kfPCRasterBase import EnKfModel
from .montecarloframework import MonteCarloModel
from .kalmanfilterframework import EnKfModel

from bmipy import Bmi
from typing import Tuple
import numpy as np

class DynamicBase(object):
    def __init__(self):
        if self.__class__ is DynamicBase:
            raise NotImplementedError

        self._d_nrTimeSteps = 0
        self.currentStep = 0
        self._d_firstTimeStep = 1
        self.inTimeStep = False
        self.inInitial = False
        self.inDynamic = False

    # def setDebug(self):
    #     msg = "Class needs to implement 'setDebug' method"
    #     raise NotImplementedError(msg)

    # def initial(self):
    #     """  """
    #     msg = "Class needs to implement 'initial' method"
    #     raise NotImplementedError(msg)

    # def dynamic(self):
    #     """  """
    #     msg = "Class needs to implement 'dynamic' method"
    #     raise NotImplementedError(msg)

    # def timeSteps(self):
    #     """  """
    #     msg = "Class needs to implement 'timeSteps' method"
    #     raise NotImplementedError(msg)

    # def nrTimeSteps(self):
    #     """  """
    #     msg = "Class needs to implement 'nrTimeSteps' method"
    #     raise NotImplementedError(msg)

    # def firstTimeStep(self):
    #     """
    #     Return the first timestep that is executed.
    #     """
    #     msg = "Class needs to implement 'firstTimeStep' method"
    #     raise NotImplementedError(msg)

    def setQuiet(self, quiet=True):
        """
        Disables the progress display of timesteps.
        """
        self.silentModelOutput = quiet

    def currentTimeStep(self):
        """  """
        msg = "Class needs to implement 'currentTimeStep' method"
        raise NotImplementedError(msg)

    def _inDynamic(self):
        return self.inDynamic

    def _inInitial(self):
        return self.inInitial

    def _setInInitial(self, value):
        assert isinstance(value, bool)
        self.inInitial = value

    def _setInDynamic(self, value):
        assert isinstance(value, bool)

        self.inDynamic = value

    def _inTimeStep(self):
        """
        Returns whether a time step is currently executing.
        """
        #if hasattr(self._userModel(), "_d_inTimeStep"):
        #  return self._userModel()._d_inTimeStep
        #else:
        #  return False
        return self.inTimeStep

    def _setInTimeStep(self, value):
        assert isinstance(value, bool)
        self.inTimeStep = value

    def _setFirstTimeStep(self, firstTimeStep):

        if not isinstance(firstTimeStep, int):
            msg = "first timestep argument (%s) of DynamicFramework must be of type int" % (type(firstTimeStep))
            raise AttributeError(msg)

        if firstTimeStep <= 0:
            msg = "first timestep argument (%s) of DynamicFramework must be > 0" % (firstTimeStep)
            raise AttributeError(msg)

        if firstTimeStep > self.nrTimeSteps():
            msg = "first timestep argument (%s) of DynamicFramework must be smaller than given last timestep (%s)" % (firstTimeStep, self.nrTimeSteps())
            raise AttributeError(msg)

        self._d_firstTimeStep = firstTimeStep

    def _setCurrentTimeStep(self, step):

        if step <= 0:
            msg = "Current timestep must be > 0"
            raise AttributeError(msg)

        if step > self.nrTimeSteps():
            msg = "Current timestep must be <= %d (nrTimeSteps)"
            raise AttributeError(msg)

        self.currentStep = step

    def _setNrTimeSteps(self, lastTimeStep):
        """
        Set the number of time steps to run.
        """
        if not isinstance(lastTimeStep, int):
            msg = "last timestep argument (%s) of DynamicFramework must be of type int" % (type(lastTimeStep))
            raise AttributeError(msg)

        if lastTimeStep <= 0:
            msg = "last timestep argument (%s) of DynamicFramework must be > 0" % (lastTimeStep)
            raise AttributeError(msg)

        self._d_nrTimeSteps = lastTimeStep

class DynamicModel(DynamicBase):

    def __init__(self):
        DynamicBase.__init__(self)
        self.silentModelOutput = False

    def initial(self):
        print("Implement 'initial' method")

    def dynamic(self):
        print("Implement 'dynamic' method")

    def _silentModelOutput(self):
        return self.silentModelOutput

    def timeSteps(self):
        return range(self.firstTimeStep(), self.nrTimeSteps() + 1)

    def nrTimeSteps(self):
        assert self._d_nrTimeSteps
        return self._d_nrTimeSteps

    def currentTimeStep(self):
        assert self.currentStep >= 0
        return self.currentStep

    def firstTimeStep(self):
        assert self._d_firstTimeStep
        return self._d_firstTimeStep

    def _setNrTimeSteps(self, timeSteps):
        DynamicBase._setNrTimeSteps(self, timeSteps)

    def _setCurrentTimeStep(self, step):
        DynamicBase._setCurrentTimeStep(self, step)

class HmDynamicBase2(DynamicModel, Bmi):

    def __init__(
            self,
            model,
            config,
            modeltime,
            domain,
            variable_list,      # TODO remove
            io = None,          # TODO
            grids = None,       # TODO
            init = None         # TODO
    ):
        DynamicModel.__init__(self)
        Bmi.__init__(self)

        self.model = model(config, modeltime, domain, init)
        # The idea is that some aspects of the configuration will be understood by the Model object, whereas other aspects will be understood by the DynamicModel object
        self.config = config
        self.variable_list = variable_list
        self.IO = io

        # TODO: `grids` should be a dictionary with the coordinates of each grid
        # GRID_DICT = {0 : ['xy'], 1 : ['x', 'y'], 2 : ['x', 'y', 'z']}  # model specific
        self.grids = grids
        # initiate the state variable object
        self.stateVar_module = stateVar(self)
        # TODO: this could be used in the case of model coupling
        try:
            self.dump_timesteps = config.DUMP['timesteps']
        except:
            self.dump_timesteps = []
        # create some flags to detect what type of simulation is being invoked
        self.is_deterministic = False
        self.apply_kalman_filter = False
        self.apply_particle_filter = False
        # # This section copied from pcraster.DynamicBase:
        # # if self.__class__ is DynamicBase:
        # #     raise NotImplementedError
        # self._d_nrTimeSteps = 0
        # self.currentStep = 0
        # self._d_firstTimeStep = 1
        # self.inTimeStep = False
        # self.inInitial = False
        # self.inDynamic = False

    # Methods defined in HmDynamicBase:
    def initiate_reporting(self, num_samples=1):
        if self.config.REPORTING['report'] == True:
            self.reporting = Reporting(self.model, self.variable_list, num_samples)
        else:
            self.reporting = DummyReporting()

    def currentSampleNumber(self):
        return 1

    def initialize(self):
        self.model.currentSampleNumber = self.currentSampleNumber()        
        self.model.time.reset()
        self.model.initial()
        self.stateVar_module.initial()
        self.reporting.initial(self.currentSampleNumber())

    def update(self):
        self.model.time.update(self.currentTimeStep())
        self.model.dynamic()
        self.stateVar_module.dynamic()
        self.reporting.dynamic(self.currentSampleNumber())

    # TODO eventually remove; this is a temporary solution
    def initial(self):
        self.initialize()

    def dynamic(self):
        self.update()

    # Methods specified by Bmi class:
    # TODO check that docstrings are correctly inherited
    def update_until(self, time: float) -> None:
        # make sure that requested end time is after now
        if time < self.current_time.timestamp():
            logging.error('`time` is prior to current model time. Please choose a new `time` and try again.')
            return
        # perform timesteps until time
        logging.info(f'Beginning simulation for {datetime.fromtimestamp(self.get_current_time()).date().isoformat()} through {datetime.fromtimestamp(time).date().isoformat()}...')
        t = timer()
        while self.get_current_time() < time:
            # if it's a new month, log a message
            current_datetime = datetime.fromtimestamp(self.get_current_time())
            if current_datetime.day == 1 and current_datetime.hour == 0:
                logging.info(f'Current model time is {current_datetime.isoformat(" ")}...')
            # advance one timestep
            self.update()
        logging.info(f'Simulation completed in {pretty_timer(timer() - t)}.')
        
    def finalize(self) -> None:
        # """Perform tear-down tasks for the model.
        pass

    def get_component_name(self) -> str:
        # """Name of the component.
        return f'hm ({self.git_hash})'  # TODO: git_hash

    def get_input_item_count(self) -> int:
        # """Count of a model's input variables.
        return len(self.IO.inputs)

    def get_output_item_count(self) -> int:
        # """Count of a model's output variables.
        return len(self.IO.outputs)

    def get_input_var_names(self) -> Tuple[str]:
        # """List of a model's input variables.
        # Input variable names must be CSDMS Standard Names, also known
        # as *long variable names*.
        # Standard Names enable the CSDMS framework to determine whether
        # an input variable in one model is equivalent to, or compatible
        # with, an output variable in another model. This allows the
        # framework to automatically connect components.
        return tuple(str(var.standard_name) for var in self.IO.inputs)

    def get_output_var_names(self) -> Tuple[str]:
        # """List of a model's output variables.
        # Output variable names must be CSDMS Standard Names, also known
        # as *long variable names*.
        return tuple(str(var.standard_name for var in self.IO.outputs))

    def get_var_grid(self, name: str) -> int:
        # """Get grid identifier for the given variable.
        return next((var.grid for var in self.IO.inputs_outputs if var.standard_name == name), None)

    def get_var_type(self, name: str) -> str:
        # """Get data type of the given variable.
        return next((var.variable_type for var in self.IO.inputs_outputs if var.standard_name == name), None)
        
    def get_var_units(self, name: str) -> str:
        # """Get units of the given variable.
        # CSDMS uses the `UDUNITS`_ standard from Unidata.
        return next((var.units for var in self.IO.inputs_outputs if var.standard_name == name), None)

    def get_var_itemsize(self, name: str) -> int:
        # """Get memory use for each array element in bytes.
        return next((var.variable_item_size for var in self.IO.inputs_outputs if var.standard_name == name), None)

    def get_var_nbytes(self, name: str) -> int:
        # """Get size, in bytes, of the given variable.
        item_size = self.get_var_itemsize(name)
        grid = self.get_var_grid(name)
        return item_size * self.get_grid_size(grid)

    def get_var_location(self, name: str) -> str:
        # """Get the grid element type that the a given variable is defined on.
        # TODO: allow the user to specify this somehow?
        return 'node'

    def get_current_time(self) -> float:
        # """Current time of the model.
        return self.model.time.current_time

    def get_start_time(self) -> float:
        # """Start time of the model.
        return self.model.time.start_time  # return pd.Timestamp?

    def get_end_time(self) -> float:
        # """End time of the model.
        return self.model.time.end_time

    def get_time_units(self) -> str:
        # """Time units of the model.
        return 's'              # TODO: check this
    
    def get_time_step(self) -> float:
        # """Current time step of the model.
        return self.model.time.time_step  # TODO: check [also, make sure in seconds]

    def get_value(self, name: str, dest: np.ndarray) -> np.ndarray:
        # """Get a copy of values of the given variable.
        var = next((var for var in self.IO.inputs_outputs if var.standard_name == name), None)
        if var in None:
            return 1
        # NB variable_class in mosartwmpy is 'state', implying that the model should maintain a 'state' attribute which contains a copy of all input/output/state variables
        dest[:] = self[var.variable_class][var.variable] # TODO: check this is an ndarray
        return 0

    def get_value_ptr(self, name: str) -> np.ndarray:
        # """Get a reference to values of the given variable.
        var = next((var for var in self.IO.inputs_outputs if var.standard_name == name), None)
        if var in None:
            raise IOError(f'Variable {name} not found in model input/output definition')
        return self[var.variable_class][var.variable]  
        
    def get_value_at_indices(
        self, name: str, dest: np.ndarray, inds: np.ndarray
    ) -> np.ndarray:
        # """Get values at particular indices.
        var = next((var for var in self.IO.inputs_outputs if var.standard_name == name), None)
        if var in None:
            return 1
        dest[:] = self[var.variable_class][var.variable][inds]
        return 0

    def set_value(self, name: str, src: np.ndarray) -> None:
        # """Specify a new value for a model variable.
        var = next((var for var in self.IO.inputs_outputs if var.standard_name == name), None)
        if var in None:
            return 1
        self[var.variable_class][var.variable][:] = src
        return 0

    def set_value_at_indices(
        self, name: str, inds: np.ndarray, src: np.ndarray
    ) -> None:
        # """Specify a new value for a model variable at particular indices.
        var = next((var for var in self.IO.inputs_outputs if var.standard_name == name), None)
        if var in None:
            return 1
        self[var.variable_class][var.variable][inds] = src
        return 0

    # Grid information
    def get_grid_rank(self, grid: int) -> int:
        # """Get number of dimensions of the computational grid.
        return len(self.grids[grid])

    def get_grid_size(self, grid: int) -> int:
        # """Get the total number of elements in the computational grid.
        coords = self.grids[grid]
        grid_size = np.prod([len(elem) for elem in coords])
        return grid_size

    def get_grid_type(self, grid: int) -> str:
        # """Get the grid type as a string.
        # Other grid types available (https://bmi-spec.readthedocs.io/en/latest/):
        # scalar
        # points
        # vector
        # unstructured
        # structured_quadrilateral
        # rectilinear
        # uniform_rectilinear
        return 'uniform_rectilinear'

    def get_grid_shape(self, grid: int, shape: np.ndarray) -> np.ndarray:
        # """Get dimensions of the computational grid.
        # This function is used for `uniform_rectilinear`,
        # `rectilinear`, `structured_quadrilateral` grids.
        coords = self.grids[grid]
        shape[:] = [len(elem) for elem in coords]
        return shape

    def get_grid_spacing(self, grid: int, spacing: np.ndarray) -> np.ndarray:
        # """Get distance between nodes of the computational grid.
        # This function is used for `uniform_rectilinear` grids.
        spacing[0] = abs(self.model.domain.y[1] - self.model.domain.y[0])
        spacing[1] = abs(self.model.domain.x[1] - self.model.domain.x[0])
        return spacing

    def get_grid_origin(self, grid: int, origin: np.ndarray) -> np.ndarray:
        # """Get coordinates for the lower-left corner of the computational grid.
        # This function is used for `uniform_rectilinear` grids.
        origin[0] = np.min(self.model.domain.y)
        origin[1] = np.min(self.model.domain.x)
        return origin

    # Non-uniform rectilinear, curvilinear
    def get_grid_x(self, grid: int, x: np.ndarray) -> np.ndarray:
        # """Get coordinates of grid nodes in the x direction.
        # This function is used for `rectilinear`, `structured_quadrilateral`,
        # and `unstructured` grids.
        x[:] = self.model.domain.x

    def get_grid_y(self, grid: int, y: np.ndarray) -> np.ndarray:
        # """Get coordinates of grid nodes in the y direction.
        # This function is used for `rectilinear`, `structured_quadrilateral`,
        # and `unstructured` grids.
        y[:] = self.model.domain.y

    def get_grid_z(self, grid: int, z: np.ndarray) -> np.ndarray:
        # """Get coordinates of grid nodes in the z direction.
        coords = self.grids[grid]
        if not 'z' in coords.keys():
            return None
        z[:] = self.model.domain.z
        return z

    def get_grid_node_count(self, grid: int) -> int:
        # """Get the number of nodes in the grid.
        # This function is used for all grid types.
        return self.model.domain.nxy

    def get_grid_edge_count(self, grid: int) -> int:
        # """Get the number of edges in the grid.
        # This function is used for `unstructured` grids.
        raise NotImplementedError

    def get_grid_face_count(self, grid: int) -> int:
        # """Get the number of faces in the grid.
        # This function is used for `unstructured` grids.
        raise NotImplementedError

    def get_grid_edge_nodes(self, grid: int, edge_nodes: np.ndarray) -> np.ndarray:
        # """Get the edge-node connectivity.
        # This function is used for `unstructured` grids.
        raise NotImplementedError

    def get_grid_face_edges(self, grid: int, face_edges: np.ndarray) -> np.ndarray:
        # """Get the face-edge connectivity.
        # This function is used for `unstructured` grids.
        raise NotImplementedError

    def get_grid_face_nodes(self, grid: int, face_nodes: np.ndarray) -> np.ndarray:
        # """Get the face-node connectivity.
        # This function is used for `unstructured` grids.
        raise NotImplementedError

    def get_grid_nodes_per_face(
        self, grid: int, nodes_per_face: np.ndarray
    ) -> np.ndarray:
        # """Get the number of nodes for each face.
        # This function is used for `unstructured` grids.
        raise NotImplementedError
    
class HmDynamicBase(DynamicModel):
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
        # # There is probably a good reason for including this here:
        # self.model.initial()
        self.config = config
        self.variable_list = variable_list
        # initiate the state variable object
        self.stateVar_module = stateVar(self)
        # TODO: this could be used in the case of model coupling
        try:
            self.dump_timesteps = config.DUMP['timesteps']
        except:
            self.dump_timesteps = []

        # create some flags to detect what type of simulation is being invoked
        self.is_deterministic = False
        self.apply_kalman_filter = False
        self.apply_particle_filter = False
        
    def initiate_reporting(self, num_samples=1):
        if self.config.REPORTING['report'] == True:
            self.reporting = Reporting(self.model, self.variable_list, num_samples)
        else:
            self.reporting = DummyReporting()

    def currentSampleNumber(self):
        return 1

    # initial() and dynamic() methods are run within the MC loop
    def initial(self):
        self.model.currentSampleNumber = self.currentSampleNumber()        
        self.model.time.reset()
        self.model.initial()
        self.stateVar_module.initial()
        self.reporting.initial(self.currentSampleNumber())

    def dynamic(self):
        self.model.time.update(self.currentTimeStep())
        self.model.dynamic()
        self.stateVar_module.dynamic()
        self.reporting.dynamic(self.currentSampleNumber())
            
class HmDynamicModel(HmDynamicBase2):
    def __init__(
        self,
        model,
        config,
        modeltime,
        domain,
        variable_list,
        init = None
    ):
        HmDynamicBase.__init__(self, model, config, modeltime, domain, variable_list, init)
        self.initiate_reporting()
        self.is_deterministic = True
        
class HmMonteCarloModel(HmDynamicBase, MonteCarloModel):
    def __init__(
        self,
        model,
        config,
        modeltime,
        domain,
        variable_list,
        init = None
    ):
        HmDynamicBase.__init__(self, model, config, modeltime, domain, variable_list, init)
        MonteCarloModel.__init__(self)

    def currentSampleNumber(self):
        return MonteCarloModel.currentSampleNumber(self)
        
    def premcloop(self):
        self.initiate_reporting(self.nrSamples())

    # NOT SURE IF NEEDED:
    # def initial(self):
    #     self.reporting.initial(self.currentSampleNumber())

    def postmcloop(self):
        pass
        
class HmEnKfModel(HmMonteCarloModel, EnKfModel):
    def __init__(
            self,
            model,
            config,
            modeltime,
            domain,
            variable_list,
            init = None
    ):
        HmMonteCarloModel.__init__(self, model, config, modeltime, domain, variable_list, init)
        EnKfModel.__init__(self)
        filter_timesteps = self.config.KALMAN_FILTER['filter_timesteps']
        combined_dump_timesteps = sorted(filter_timesteps + self.dump_timesteps)
        self.dump_timesteps = combined_dump_timesteps        
        # EXPERIMENTAL:
        self.model.filter_timesteps = filter_timesteps

    def setState(self):
        pass

    def setObservations(self):
        pass

    def resume(self):
        self.stateVar_module.resume()

class HmDynamicCoupledModel(HmMonteCarloModel):
    def __init__(
            self,
            model,
            config,
            modeltime,
            domain,
            variable_list,
            init = None
    ):
        HmMonteCarloModel.__init__(self, model, config, modeltime, domain, variable_list, init)
        # EnKfModel.__init__(self)
        # filter_timesteps = self.config.KALMAN_FILTER['filter_timesteps']
        # combined_dump_timesteps = sorted(filter_timesteps + self.dump_timesteps)
        # self.dump_timesteps = combined_dump_timesteps

    def setState(self):
        pass

    def setObservations(self):
        pass

    def resume(self):
        self.stateVar_module.resume()
        
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
#         DynamicModel.__init__(self)
#         self.modeltime = modeltime
#         self.model = model(
#             config,
#             modeltime,
#             domain,
#             init
#         )
#         self.model.initial()
#         if config.REPORTING['report'] == True:
#             self.reporting = Reporting(self.model, variable_list)
#         else:
#             self.reporting = DummyReporting()

#     def initial(self):
#         self.reporting.initial()

#     def dynamic(self):
#         self.model.time.update(self.currentTimeStep())
#         self.model.dynamic()
#         self.reporting.dynamic()
    
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
