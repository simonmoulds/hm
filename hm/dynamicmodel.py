#!/usr/bin/env python
# -*- coding: utf-8 -*-

from .reporting import Reporting, DummyReporting
from .stateVar import stateVar
from .pcraster.dynamicPCRasterBase import DynamicModel
from .pcraster.mcPCRasterBase import MonteCarloModel
from .pcraster.kfPCRasterBase import EnKfModel

from bmipy import bmi

# class DynamicBase(object):
#   """
#   Base class for DynamicModel framework class.

#   .. todo::

#     KDJ: What is the reason this class exists? Can it be merged with
#     DynamicModel?
#   """

#   def __init__(self):
#     if self.__class__ is DynamicBase:
#       raise NotImplementedError

#     self._d_nrTimeSteps = 0
#     self.currentStep = 0
#     self._d_firstTimeStep = 1
#     self.inTimeStep = False
#     self.inInitial = False
#     self.inDynamic = False

#   def setDebug(self):
#     msg = "Class needs to implement 'setDebug' method"
#     raise NotImplementedError(msg)

#   def initial(self):
#     """  """
#     msg = "Class needs to implement 'initial' method"
#     raise NotImplementedError(msg)

#   def dynamic(self):
#     """  """
#     msg = "Class needs to implement 'dynamic' method"
#     raise NotImplementedError(msg)

#   def timeSteps(self):
#     """  """
#     msg = "Class needs to implement 'timeSteps' method"
#     raise NotImplementedError(msg)

#   def nrTimeSteps(self):
#     """  """
#     msg = "Class needs to implement 'nrTimeSteps' method"
#     raise NotImplementedError(msg)

#   def firstTimeStep(self):
#     """
#     Return the first timestep that is executed.
#     """
#     msg = "Class needs to implement 'firstTimeStep' method"
#     raise NotImplementedError(msg)

#   def setQuiet(self,
#     quiet=True):
#     """
#     Disables the progress display of timesteps.
#     """
#     self.silentModelOutput = quiet

#   def currentTimeStep(self):
#     """  """
#     msg = "Class needs to implement 'currentTimeStep' method"
#     raise NotImplementedError(msg)

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

#   def _setCurrentTimeStep(self, step):

#     if step <= 0:
#       msg = "Current timestep must be > 0"
#       raise AttributeError(msg)

#     if step > self.nrTimeSteps():
#       msg = "Current timestep must be <= %d (nrTimeSteps)"
#       raise AttributeError(msg)

#     self.currentStep = step

#   def _setNrTimeSteps(self, lastTimeStep):
#     """
#     Set the number of time steps to run.
#     """
#     if not isinstance(lastTimeStep, int):
#       msg = "last timestep argument (%s) of DynamicFramework must be of type int" % (type(lastTimeStep))
#       raise AttributeError(msg)

#     if lastTimeStep <= 0:
#       msg = "last timestep argument (%s) of DynamicFramework must be > 0" % (lastTimeStep)
#       raise AttributeError(msg)

#     self._d_nrTimeSteps = lastTimeStep

# class DynamicModel(dynamicBase.DynamicBase):

#   def __init__(self):
#     dynamicBase.DynamicBase.__init__(self)
#     self.silentModelOutput = False

#   def initial(self):
#     print("Implement 'initial' method")

#   def dynamic(self):
#     print("Implement 'dynamic' method")

#   def _silentModelOutput(self):
#     return self.silentModelOutput

#   def timeSteps(self):
#     """
#     Return a list of time steps configured
#     """
#     return range(self.firstTimeStep(), self.nrTimeSteps() + 1)

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

#   def report(self,
#     variable,
#     name):
#     """
#     Storing map data to disk

#     `variable`
#       Variable containing the PCRaster map data

#     `name`
#       Name used as filename. Use a filename with less than eight
#       characters and without extension. File extension for dynamic models
#       is ".map" in the initial section and the 8.3 style format name in
#       the dynamic section. File extensions will be appended automatically.
#     """
#     self._reportNew(variable, name)

#   def readmap(self, name, style=1):
#     """
#     Read map data from disk.

#     `name`
#       Name used as filename. Use filename with less than eight characters
#       and without extension. File extension for dynamic models is ".map"
#       in initial section and the 8.3 style format name in the dynamic
#       section. File extensions will be appended automatically.

#     .. todo::

#       `style` argument is not used.
#     """
#     return self._readmapNew(name)

#   def _setNrTimeSteps(self,
#     timeSteps):
#     """
#     Configure the number of time steps.

#     In addition to the setting the number of timesteps we need to pass
#     the value to the PCRaster runtime engine.
#     """
#     dynamicBase.DynamicBase._setNrTimeSteps(self, timeSteps)

#     # pcr._rte().setNrTimeSteps(timeSteps)

#   def _setCurrentTimeStep(self,
#     step):
#     """
#     Set the current time step.

#     In addition to the setting the current timestep within the framework,
#     we need to pass the value to the PCRaster runtime engine.
#     """
#     dynamicBase.DynamicBase._setCurrentTimeStep(self, step)

#     # pcr._rte().setCurrentTimeStep(step)

class HmDynamicBase2(Bmi):

    def __init__(
            self,
            model,
            config,
            modeltime,
            domain,
            variable_list,
            init = None
    ):
        self.model = model(config, modeltime, domain, init)
        # # There is probably a good reason for including this here:
        # self.model.initial()

        # The idea is that some aspects of the configuration will be understood by the Model object, whereas other aspects will be understood by the DynamicModel object
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
        # self.apply_kalman_filter = False
        # self.apply_particle_filter = False
            
    def initiate_reporting(self, num_samples=1):
        if self.config.REPORTING['report'] == True:
            self.reporting = Reporting(self.model, self.variable_list, num_samples)
        else:
            self.reporting = DummyReporting()

    # def currentSampleNumber(self):
    #     return 1
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

    # Now implement all the other methods specified by the Bmi class
    
    # def update_until(self, time: float) -> None:
    #     """Advance model state until the given time.

    #     Parameters
    #     ----------
    #     time : float
    #         A model time later than the current model time.
    #     """        
    #     ...

    # @abstractmethod
    # def finalize(self) -> None:
    #     """Perform tear-down tasks for the model.

    #     Perform all tasks that take place after exiting the model's time
    #     loop. This typically includes deallocating memory, closing files and
    #     printing reports.
    #     """
    #     ...

    # @abstractmethod
    # def get_component_name(self) -> str:
    #     """Name of the component.

    #     Returns
    #     -------
    #     str
    #         The name of the component.
    #     """
    #     ...

    # @abstractmethod
    # def get_input_item_count(self) -> int:
    #     """Count of a model's input variables.

    #     Returns
    #     -------
    #     int
    #       The number of input variables.
    #     """
    #     ...

    # @abstractmethod
    # def get_output_item_count(self) -> int:
    #     """Count of a model's output variables.

    #     Returns
    #     -------
    #     int
    #       The number of output variables.
    #     """
    #     ...

    # @abstractmethod
    # def get_input_var_names(self) -> Tuple[str]:
    #     """List of a model's input variables.

    #     Input variable names must be CSDMS Standard Names, also known
    #     as *long variable names*.

    #     Returns
    #     -------
    #     list of str
    #         The input variables for the model.

    #     Notes
    #     -----
    #     Standard Names enable the CSDMS framework to determine whether
    #     an input variable in one model is equivalent to, or compatible
    #     with, an output variable in another model. This allows the
    #     framework to automatically connect components.

    #     Standard Names do not have to be used within the model.
    #     """
    #     ...

    # @abstractmethod
    # def get_output_var_names(self) -> Tuple[str]:
    #     """List of a model's output variables.

    #     Output variable names must be CSDMS Standard Names, also known
    #     as *long variable names*.

    #     Returns
    #     -------
    #     list of str
    #         The output variables for the model.
    #     """
    #     ...

    # @abstractmethod
    # def get_var_grid(self, name: str) -> int:
    #     """Get grid identifier for the given variable.

    #     Parameters
    #     ----------
    #     name : str
    #         An input or output variable name, a CSDMS Standard Name.

    #     Returns
    #     -------
    #     int
    #       The grid identifier.
    #     """
    #     ...

    # @abstractmethod
    # def get_var_type(self, name: str) -> str:
    #     """Get data type of the given variable.

    #     Parameters
    #     ----------
    #     name : str
    #         An input or output variable name, a CSDMS Standard Name.

    #     Returns
    #     -------
    #     str
    #         The Python variable type; e.g., ``str``, ``int``, ``float``.
    #     """
    #     ...

    # @abstractmethod
    # def get_var_units(self, name: str) -> str:
    #     """Get units of the given variable.

    #     Standard unit names, in lower case, should be used, such as
    #     ``meters`` or ``seconds``. Standard abbreviations, like ``m`` for
    #     meters, are also supported. For variables with compound units,
    #     each unit name is separated by a single space, with exponents
    #     other than 1 placed immediately after the name, as in ``m s-1``
    #     for velocity, ``W m-2`` for an energy flux, or ``km2`` for an
    #     area.

    #     Parameters
    #     ----------
    #     name : str
    #         An input or output variable name, a CSDMS Standard Name.

    #     Returns
    #     -------
    #     str
    #         The variable units.

    #     Notes
    #     -----
    #     CSDMS uses the `UDUNITS`_ standard from Unidata.

    #     .. _UDUNITS: http://www.unidata.ucar.edu/software/udunits
    #     """
    #     ...

    # @abstractmethod
    # def get_var_itemsize(self, name: str) -> int:
    #     """Get memory use for each array element in bytes.

    #     Parameters
    #     ----------
    #     name : str
    #         An input or output variable name, a CSDMS Standard Name.

    #     Returns
    #     -------
    #     int
    #         Item size in bytes.
    #     """
    #     ...

    # @abstractmethod
    # def get_var_nbytes(self, name: str) -> int:
    #     """Get size, in bytes, of the given variable.

    #     Parameters
    #     ----------
    #     name : str
    #         An input or output variable name, a CSDMS Standard Name.

    #     Returns
    #     -------
    #     int
    #         The size of the variable, counted in bytes.
    #     """
    #     ...

    # @abstractmethod
    # def get_var_location(self, name: str) -> str:
    #     """Get the grid element type that the a given variable is defined on.

    #     The grid topology can be composed of *nodes*, *edges*, and *faces*.

    #     *node*
    #         A point that has a coordinate pair or triplet: the most
    #         basic element of the topology.

    #     *edge*
    #         A line or curve bounded by two *nodes*.

    #     *face*
    #         A plane or surface enclosed by a set of edges. In a 2D
    #         horizontal application one may consider the word “polygon”,
    #         but in the hierarchy of elements the word “face” is most common.

    #     Parameters
    #     ----------
    #     name : str
    #         An input or output variable name, a CSDMS Standard Name.

    #     Returns
    #     -------
    #     str
    #         The grid location on which the variable is defined. Must be one of
    #         `"node"`, `"edge"`, or `"face"`.

    #     Notes
    #     -----
    #     CSDMS uses the `ugrid conventions`_ to define unstructured grids.

    #     .. _ugrid conventions: http://ugrid-conventions.github.io/ugrid-conventions
    #     """
    #     ...

    # @abstractmethod
    # def get_current_time(self) -> float:
    #     """Current time of the model.

    #     Returns
    #     -------
    #     float
    #         The current model time.
    #     """
    #     ...

    # @abstractmethod
    # def get_start_time(self) -> float:
    #     """Start time of the model.

    #     Model times should be of type float.

    #     Returns
    #     -------
    #     float
    #         The model start time.
    #     """
    #     ...

    # @abstractmethod
    # def get_end_time(self) -> float:
    #     """End time of the model.

    #     Returns
    #     -------
    #     float
    #         The maximum model time.
    #     """
    #     ...

    # @abstractmethod
    # def get_time_units(self) -> str:
    #     """Time units of the model.

    #     Returns
    #     -------
    #     str
    #         The model time unit; e.g., `days` or `s`.

    #     Notes
    #     -----
    #     CSDMS uses the UDUNITS standard from Unidata.
    #     """
    #     ...

    # @abstractmethod
    # def get_time_step(self) -> float:
    #     """Current time step of the model.

    #     The model time step should be of type float.

    #     Returns
    #     -------
    #     float
    #         The time step used in model.
    #     """
    #     ...

    # @abstractmethod
    # def get_value(self, name: str, dest: np.ndarray) -> np.ndarray:
    #     """Get a copy of values of the given variable.

    #     This is a getter for the model, used to access the model's
    #     current state. It returns a *copy* of a model variable, with
    #     the return type, size and rank dependent on the variable.

    #     Parameters
    #     ----------
    #     name : str
    #         An input or output variable name, a CSDMS Standard Name.
    #     dest : ndarray
    #         A numpy array into which to place the values.

    #     Returns
    #     -------
    #     ndarray
    #         The same numpy array that was passed as an input buffer.
    #     """
    #     ...

    # @abstractmethod
    # def get_value_ptr(self, name: str) -> np.ndarray:
    #     """Get a reference to values of the given variable.

    #     This is a getter for the model, used to access the model's
    #     current state. It returns a reference to a model variable,
    #     with the return type, size and rank dependent on the variable.

    #     Parameters
    #     ----------
    #     name : str
    #         An input or output variable name, a CSDMS Standard Name.

    #     Returns
    #     -------
    #     array_like
    #         A reference to a model variable.
    #     """
    #     ...

    # @abstractmethod
    # def get_value_at_indices(
    #     self, name: str, dest: np.ndarray, inds: np.ndarray
    # ) -> np.ndarray:
    #     """Get values at particular indices.

    #     Parameters
    #     ----------
    #     name : str
    #         An input or output variable name, a CSDMS Standard Name.
    #     dest : ndarray
    #         A numpy array into which to place the values.
    #     inds : array_like
    #         The indices into the variable array.

    #     Returns
    #     -------
    #     array_like
    #         Value of the model variable at the given location.
    #     """
    #     ...

    # @abstractmethod
    # def set_value(self, name: str, src: np.ndarray) -> None:
    #     """Specify a new value for a model variable.

    #     This is the setter for the model, used to change the model's
    #     current state. It accepts, through *src*, a new value for a
    #     model variable, with the type, size and rank of *src*
    #     dependent on the variable.

    #     Parameters
    #     ----------
    #     name : str
    #         An input or output variable name, a CSDMS Standard Name.
    #     src : array_like
    #         The new value for the specified variable.
    #     """
    #     ...

    # @abstractmethod
    # def set_value_at_indices(
    #     self, name: str, inds: np.ndarray, src: np.ndarray
    # ) -> None:
    #     """Specify a new value for a model variable at particular indices.

    #     Parameters
    #     ----------
    #     name : str
    #         An input or output variable name, a CSDMS Standard Name.
    #     inds : array_like
    #         The indices into the variable array.
    #     src : array_like
    #         The new value for the specified variable.
    #     """
    #     ...

    # # Grid information
    # @abstractmethod
    # def get_grid_rank(self, grid: int) -> int:
    #     """Get number of dimensions of the computational grid.

    #     Parameters
    #     ----------
    #     grid : int
    #         A grid identifier.

    #     Returns
    #     -------
    #     int
    #         Rank of the grid.
    #     """
    #     ...

    # @abstractmethod
    # def get_grid_size(self, grid: int) -> int:
    #     """Get the total number of elements in the computational grid.

    #     Parameters
    #     ----------
    #     grid : int
    #         A grid identifier.

    #     Returns
    #     -------
    #     int
    #         Size of the grid.
    #     """
    #     ...

    # @abstractmethod
    # def get_grid_type(self, grid: int) -> str:
    #     """Get the grid type as a string.

    #     Parameters
    #     ----------
    #     grid : int
    #         A grid identifier.

    #     Returns
    #     -------
    #     str
    #         Type of grid as a string.
    #     """
    #     ...

    # # Uniform rectilinear
    # @abstractmethod
    # def get_grid_shape(self, grid: int, shape: np.ndarray) -> np.ndarray:
    #     """Get dimensions of the computational grid.

    #     Parameters
    #     ----------
    #     grid : int
    #         A grid identifier.
    #     shape : ndarray of int, shape *(ndim,)*
    #         A numpy array into which to place the shape of the grid.

    #     Returns
    #     -------
    #     ndarray of int
    #         The input numpy array that holds the grid's shape.
    #     """
    #     ...

    # @abstractmethod
    # def get_grid_spacing(self, grid: int, spacing: np.ndarray) -> np.ndarray:
    #     """Get distance between nodes of the computational grid.

    #     Parameters
    #     ----------
    #     grid : int
    #         A grid identifier.
    #     spacing : ndarray of float, shape *(ndim,)*
    #         A numpy array to hold the spacing between grid rows and columns.

    #     Returns
    #     -------
    #     ndarray of float
    #         The input numpy array that holds the grid's spacing.
    #     """
    #     ...

    # @abstractmethod
    # def get_grid_origin(self, grid: int, origin: np.ndarray) -> np.ndarray:
    #     """Get coordinates for the lower-left corner of the computational grid.

    #     Parameters
    #     ----------
    #     grid : int
    #         A grid identifier.
    #     origin : ndarray of float, shape *(ndim,)*
    #         A numpy array to hold the coordinates of the lower-left corner of
    #         the grid.

    #     Returns
    #     -------
    #     ndarray of float
    #         The input numpy array that holds the coordinates of the grid's
    #         lower-left corner.
    #     """
    #     ...

    # # Non-uniform rectilinear, curvilinear
    # @abstractmethod
    # def get_grid_x(self, grid: int, x: np.ndarray) -> np.ndarray:
    #     """Get coordinates of grid nodes in the x direction.

    #     Parameters
    #     ----------
    #     grid : int
    #         A grid identifier.
    #     x : ndarray of float, shape *(nrows,)*
    #         A numpy array to hold the x-coordinates of the grid node columns.

    #     Returns
    #     -------
    #     ndarray of float
    #         The input numpy array that holds the grid's column x-coordinates.
    #     """
    #     ...

    # @abstractmethod
    # def get_grid_y(self, grid: int, y: np.ndarray) -> np.ndarray:
    #     """Get coordinates of grid nodes in the y direction.

    #     Parameters
    #     ----------
    #     grid : int
    #         A grid identifier.
    #     y : ndarray of float, shape *(ncols,)*
    #         A numpy array to hold the y-coordinates of the grid node rows.

    #     Returns
    #     -------
    #     ndarray of float
    #         The input numpy array that holds the grid's row y-coordinates.
    #     """
    #     ...

    # @abstractmethod
    # def get_grid_z(self, grid: int, z: np.ndarray) -> np.ndarray:
    #     """Get coordinates of grid nodes in the z direction.

    #     Parameters
    #     ----------
    #     grid : int
    #         A grid identifier.
    #     z : ndarray of float, shape *(nlayers,)*
    #         A numpy array to hold the z-coordinates of the grid nodes layers.

    #     Returns
    #     -------
    #     ndarray of float
    #         The input numpy array that holds the grid's layer z-coordinates.
    #     """
    #     ...

    # @abstractmethod
    # def get_grid_node_count(self, grid: int) -> int:
    #     """Get the number of nodes in the grid.

    #     Parameters
    #     ----------
    #     grid : int
    #         A grid identifier.

    #     Returns
    #     -------
    #     int
    #         The total number of grid nodes.
    #     """
    #     ...

    # @abstractmethod
    # def get_grid_edge_count(self, grid: int) -> int:
    #     """Get the number of edges in the grid.

    #     Parameters
    #     ----------
    #     grid : int
    #         A grid identifier.

    #     Returns
    #     -------
    #     int
    #         The total number of grid edges.
    #     """
    #     ...

    # @abstractmethod
    # def get_grid_face_count(self, grid: int) -> int:
    #     """Get the number of faces in the grid.

    #     Parameters
    #     ----------
    #     grid : int
    #         A grid identifier.

    #     Returns
    #     -------
    #     int
    #         The total number of grid faces.
    #     """
    #     ...

    # @abstractmethod
    # def get_grid_edge_nodes(self, grid: int, edge_nodes: np.ndarray) -> np.ndarray:
    #     """Get the edge-node connectivity.

    #     Parameters
    #     ----------
    #     grid : int
    #         A grid identifier.
    #     edge_nodes : ndarray of int, shape *(2 x nnodes,)*
    #         A numpy array to place the edge-node connectivity. For each edge,
    #         connectivity is given as node at edge tail, followed by node at
    #         edge head.

    #     Returns
    #     -------
    #     ndarray of int
    #         The input numpy array that holds the edge-node connectivity.
    #     """
    #     ...

    # @abstractmethod
    # def get_grid_face_edges(self, grid: int, face_edges: np.ndarray) -> np.ndarray:
    #     """Get the face-edge connectivity.

    #     Parameters
    #     ----------
    #     grid : int
    #         A grid identifier.
    #     face_edges : ndarray of int
    #         A numpy array to place the face-edge connectivity.

    #     Returns
    #     -------
    #     ndarray of int
    #         The input numpy array that holds the face-edge connectivity.
    #     """
    #     ...

    # @abstractmethod
    # def get_grid_face_nodes(self, grid: int, face_nodes: np.ndarray) -> np.ndarray:
    #     """Get the face-node connectivity.

    #     Parameters
    #     ----------
    #     grid : int
    #         A grid identifier.
    #     face_nodes : ndarray of int
    #         A numpy array to place the face-node connectivity. For each face,
    #         the nodes (listed in a counter-clockwise direction) that form the
    #         boundary of the face.

    #     Returns
    #     -------
    #     ndarray of int
    #         The input numpy array that holds the face-node connectivity.
    #     """
    #     ...

    # @abstractmethod
    # def get_grid_nodes_per_face(
    #     self, grid: int, nodes_per_face: np.ndarray
    # ) -> np.ndarray:
    #     """Get the number of nodes for each face.

    #     Parameters
    #     ----------
    #     grid : int
    #         A grid identifier.
    #     nodes_per_face : ndarray of int, shape *(nfaces,)*
    #         A numpy array to place the number of nodes per face.

    #     Returns
    #     -------
    #     ndarray of int
    #         The input numpy array that holds the number of nodes per face.
    #     """
    #     ...
    
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
            
class HmDynamicModel(HmDynamicBase):
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
