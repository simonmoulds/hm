# !/usr/bin/env python
# -*- coding: utf-8 -*-

# from .pcraster.dynamicFramework import DynamicFramework
from .pcraster.mcFramework import MonteCarloFramework
from .pcraster.kalmanFilterFramework import EnsKalmanFilterFramework

import atexit
import os
import re
import sys
import weakref
# try:
#   from . import shellscript
# except ImportError as error:
#   print("PCRaster modelling framework error: {}".format(error))
#   raise SystemExit

class FrameworkError(Exception):
  def __init__(self,
    msg):
    self._msg = msg

  def __str__(self):
    return self._msg

# class WeakCallback (object):
#     """A Weak Callback object that will keep a reference to
#     the connecting object with weakref semantics.

#     This allows object A to pass a callback method to object S,
#     without object S keeping A alive.
#     """
#     def __init__(self, mcallback):
#         """Create a new Weak Callback calling the method @mcallback"""
#         if sys.version_info.major == 2:
#             obj = mcallback.im_self
#             attr = mcallback.im_func.__name__
#         else:
#             obj = mcallback.__self__
#             attr = mcallback.__func__.__name__

#         self.wref = weakref.ref(obj, self.object_deleted)
#         self.callback_attr = attr

#     def __call__(self, *args, **kwargs):
#         obj = self.wref()
#         if obj:
#             attr = getattr(obj, self.callback_attr)
#             result = attr(*args, **kwargs)
#         else:
#             result = self.default_callback(*args, **kwargs)
#         return result

#     def default_callback(self, *args, **kwargs):
#         """Called instead of callback when expired"""
#         assert False
#         # pass

#     def object_deleted(self, wref):
#         """Called when callback expires"""
#         pass

# @atexit.register
# def _atExit():
#   print

# def generateNameT(
#   name,
#   time):
#   """
#   Return a filename based on the name and time step passed in.

#   The resulting name obeys the 8.3 DOS style format. The time step
#   will be added to the end of the filename and be prepended by 0's if
#   needed.

#   The time step normally ranges from [1, nrTimeSteps].
#   The length of the name should be max 8 characters to leave room for
#   the time step.

#   The name passed in may contain a directory name.

#   See also: generateNameS(), generateNameST()
#   """
#   head, tail = os.path.split(name)

#   if re.search("\.", tail):
#     msg = "File extension given in '" + name + "' not allowed"
#     raise FrameworkError(msg)

#   if len(tail) == 0:
#     msg = "No filename specified"
#     raise FrameworkError(msg)

#   if len(tail) > 8:
#     msg = "Filename '" + name + "' must be shorter than 8 characters"
#     raise FrameworkError(msg)

#   if time < 0:
#     msg = "Timestep must be larger than 0"
#     raise FrameworkError(msg)

#   nr = "%d" % (time)

#   space = 11 - (len(tail) + len(nr))
#   assert space >= 0

#   result = "%s%s%s" % (tail, space * "0", nr)
#   result = "%s.%s" % (result[:8], result[8:])
#   assert len(result) == 12

#   return os.path.join(head, result)



# def generateNameS(name, sample):
#   """Return a filename based on the name and sample number passed in.

#   The resulting name contains a directory and a filename part. The
#   sample number is used as the directory name and the name is used as
#   the filename.

#   The sample number normally ranges from [1, nrSamples].

#   See also: generateNameT(), generateNameST()
#   """
#   return os.path.join("%d" % (sample), name)



# def generateNameST(name, sample, timestep):
#   """
#   Return a filename based on the name, sample number and time step.

#   See also: generateNameT(), generateNameS()
#   """
#   return generateNameS(generateNameT(name, timestep), sample)

# class FrameworkBase(shellscript.ShellScript):
class FrameworkBase(object):
    """
    Base class for frameworks.

    Basically contains things for logging...
    """
    # output relevant attributes
    _d_quiet = False
    _d_trace = False
    _d_debug = False
    _d_indentLevel = 0
    _d_inScript = False
    # _d_mapExtension = ".map"
    def __init__(self):
        # shellscript.ShellScript.__init__(self, ["frameworkBase.py"])
        self._d_silentModelOutput = False
        self._d_silentFrameworkOutput = True
        self._d_quietProgressDots = False
        self._d_quietProgressSampleNr = False

    def _inUpdateWeight(self):
        if not hasattr(self._userModel(), "_d_inUpdateWeight"):
            return False
        return self._userModel()._d_inUpdateWeight

    def _inResume(self):
        if not hasattr(self._userModel(), "_d_inResume"):
            return False
        return self._userModel()._d_inResume

    def __del__(self):
        self._atEndOfScript()

    def _silentModelOutput(self):
        return self._d_silentModelOutput

    # def _addMethodToClass(self, func):
    #     # NO! This will create a circular reference between the user model and the
    #     # framework class. Both will never be deleted because the reference counts
    #     # will never drop to zero.
    #     # setattr(self._userModel(), func.__name__, func)

    #     # Use a weak reference to the framework class. This assumes that the
    #     # framework class will be alive for as long as the user model is used.
    #     call_back = WeakCallback(func)
    #     setattr(self._userModel(), func.__name__, call_back)

    # def _addAttributeToClass(self, attr, value):
    #     setattr(self._userModel(), attr, value)

    def _atStartOfTimeStep(self, step):
        self._userModel()._setInTimeStep(True)
        if not self._quiet():
            if not self._trace():
                msg = u"."
            else:
                msg = u"%s<time step=\"%s\">\n" % (self._indentLevel(), step)
                sys.stdout.write(msg)
                sys.stdout.flush()

    def _timeStepFinished(self):
        self._userModel()._setInTimeStep(False)
        if not self._quiet():
            if self._trace():
                self.showMessage("%s</time>" % (self._indentLevel()))

    def _atStartOfScript(self):
        if not self._d_inScript:
            self._userModel()._d_inScript = True
            if not self._quiet():
                if self._trace():
                    self.showMessage("<script>")

    def _atEndOfScript(self):
        if self._d_inScript:
            self._d_inScript = False
            if not self._quiet():
                if not self._trace():
                    msg = u"\n"
                else:
                    msg = u"</script>\n"
                # showMessage does not work due to encode throw
                sys.stdout.write(msg)
                sys.stdout.flush()

    def _incrementIndentLevel(self):
        FrameworkBase._d_indentLevel += 1

    def _decrementIndentLevel(self):
        assert FrameworkBase._d_indentLevel > 0
        FrameworkBase._d_indentLevel -= 1

    def _scriptFinished(self):
        self._atEndOfScript()

    def _trace(self):
        return FrameworkBase._d_trace

    def _debug(self):
        return FrameworkBase._d_debug

    def _indentLevel(self):
        return FrameworkBase._d_indentLevel * "  "

    def _traceIn(self, functionName):
        if not self._quiet():
            if self._trace():
                self.showMessage("%s<%s>" % (self._indentLevel(), functionName))

    def _traceOut(self, functionName):
        if not self._quiet():
            if self._trace():
                self.showMessage("%s</%s>" % (self._indentLevel(), functionName))

    def _quiet(self):
        return self._d_quietProgressDots

    def _quietProgressSampleNr(self):
        """Return state of sample number display.
        .. todo::
        This method assumes a Monte Carlo style framework specialization.
        We should think of a more general verbosity level member to which
        frameworks can respond.
        """
        return self._d_quietProgressSampleNr

    def setQuiet(self, quiet):
        """
        Enable/disable all framework output to stdout.
        `quiet`
        True/False. Default is set to False
        """
        FrameworkBase._d_quiet = quiet

    def setTrace(self, trace):
        """
        Trace framework output to stdout.
        `trace`
        True/False. Default is set to False.

        If tracing is enabled the user will get a detailed framework output
        in an XML style.
        """
        FrameworkBase._d_trace = trace

    def setDebug(self, debug):
        FrameworkBase._d_debug = debug

    def _atStartOfSample(self, nr):
        self._userModel()._d_inSample = True
        if not self._quietProgressSampleNr():
            if not self._trace():
                msg = u"%d " % (nr)
            else:
                msg = u"%s<sample nr=\"%s\">\n" % (self._indentLevel(), nr)
            # no showMessage here, \n not desired in non-trace "..." timestep output
            sys.stdout.write(msg)
            sys.stdout.flush()

    def _sampleFinished(self):
        self._userModel()._d_inSample = False
        if not self._quiet():
            #if not self._trace():
            #msg = "]"
            #else:
            if self._trace():
                msg = "%s</sample>" % (self._indentLevel())
                self.showMessage(msg)

    def _atStartOfFilterPeriod(self, nr):
        self._userModel()._d_inFilterPeriod = True
        if not self._d_model._quiet():
            if not self._d_model._trace():
                msg = "\nPeriod %d" % (nr + 1)
            else:
                msg = "%s<period nr=\"%s\">" % (self._indentLevel(), nr + 1)
            self.showMessage(msg)

    def _atEndOfFilterPeriod(self):
        self._userModel()._d_inFilterPeriod = False
        if not self._d_model._quiet():
            if self._d_model._trace():
                msg = "%s</period>" % (self._indentLevel())
                self.showMessage(msg)

    def _runInitial(self):
        self._userModel()._setInInitial(True)
        if(hasattr(self._userModel(), 'initial')):
            self._incrementIndentLevel()
            self._traceIn("initial")
            self._userModel().initial()
            self._traceOut("initial")
            self._decrementIndentLevel()
        self._userModel()._setInInitial(False)

    def _runDynamic(self):
        self._userModel()._setInDynamic(True)
        step = self._userModel().firstTimeStep()
        while step <= self._userModel().nrTimeSteps():
            self._incrementIndentLevel()
            self._atStartOfTimeStep(step)
            self._userModel()._setCurrentTimeStep(step)
            if hasattr(self._userModel(), 'dynamic'):
                self._incrementIndentLevel()
                self._traceIn("dynamic")
                self._userModel().dynamic()
                self._traceOut("dynamic")
                self._decrementIndentLevel()
            self._timeStepFinished()
            self._decrementIndentLevel()
            step += 1
        self._userModel()._setInDynamic(False)

    def _runSuspend(self):
        if(hasattr(self._userModel(), 'suspend')):
            self._incrementIndentLevel()
            self._traceIn("suspend")
            self._userModel().suspend()
            self._traceOut("suspend")
            self._decrementIndentLevel()

    def _runResume(self):
        self._userModel()._d_inResume = True
        if(hasattr(self._userModel(), 'resume')):
            self._incrementIndentLevel()
            self._traceIn("resume")
            self._userModel().resume()
            self._traceOut("resume")
            self._decrementIndentLevel()
        self._userModel()._d_inResume = False

    def _runPremcloop(self):
        self._userModel()._d_inPremc = True
        if hasattr(self._userModel(), 'premcloop') :
            self._incrementIndentLevel()
            self._traceIn("premcloop")
            self._userModel().premcloop()
            self._traceOut("premcloop")
            self._decrementIndentLevel()
        self._userModel()._d_inPremc = False

    def _runPostmcloop(self):
        self._userModel()._d_inPostmc = True
        if hasattr(self._userModel(), 'postmcloop'):
            self._incrementIndentLevel()
            self._traceIn("postmcloop")
            self._userModel().postmcloop()
            self._traceOut("postmcloop")
            self._decrementIndentLevel()
        self._userModel()._d_inPostmc = False

    # def _report(self, variable, name):
    #     """
    #     Report map data to disk.

    #     .. todo::

    #       Uses PCRaster package which isn't imported.
    #     """
    #     pass
    #     # pcraster.report(variable, name)

    # def _generateName(self, name):
    #   """
    #   Return a filename based on the name and the current sample number
    #   and time step.

    #   The filename created will depend on whether this function is called
    #   from within a sample and/or a time step. Pseudo code:

    #   if sample and time step:
    #     sample/name.timestep
    #   else if sample:
    #     sample/name
    #   else if time step:
    #     name.timestep
    #   else:
    #     name

    #   If this function is not called from within a time step and name does not
    #   have an extension, the default extension '.map' is added to the name.

    #   See also: generateNameS(), generateNameT(), generateNameST()
    #   """
    #   if self._inSample() and self._inTimeStep():
    #     name = self._generateNameST(name, self.currentSampleNumber(),
    #       self.currentTimeStep())
    #   elif self.inSample():
    #     name = self._generateNameS(name, self.currentSampleNumber())
    #   elif self.inTimeStep():
    #     name = self._generateNameT(name, self.currentTimeStep())

    #   if not self.inTimeStep() and len(os.path.splitext(name)[1]) == 0:
    #     name = name + ".map"

    #   return name

    # def _generateNameT(self, name, time):
    #   return generateNameT(name, time)

    # def _generateNameS(self, name, sample):
    #   return generateNameS(name, sample)

    # def _generateNameST(self, name, sample, time):
    #   """
    #   Return a filename based on the name, sample number and time step.

    #   See also: generateNameT(), generateNameS()
    #   """
    #   return self._generateNameS(self._generateNameT(name, time), sample)

    # def generateNameS(self, name, sample):
    #   return generateNameS(name, sample)

    # def _reportNew(self,
    #   variable,
    #   name,
    #   style=1):
    #   """

    #   .. todo::

    #     `style` argument is not used.
    #   """
    #   head, tail = os.path.split(name)

    #   if re.search("\.", tail):
    #     msg = "File extension given in '" + name + "' not allowed, provide filename without extension"
    #     raise FrameworkError(msg)

    #   directoryPrefix = ""
    #   nameSuffix = ".map"
    #   newName = ""

    #   if hasattr(self._userModel(), "_inStochastic"):
    #     if self._userModel()._inStochastic():
    #       if self._userModel()._inPremc():
    #         newName = name + nameSuffix
    #       elif self._userModel()._inPostmc():
    #         newName = name + nameSuffix
    #       else:
    #         directoryPrefix = str(self._userModel().currentSampleNumber())

    #   if self._userModel()._inInitial():
    #     newName = name + nameSuffix

    #   if hasattr(self._userModel(), "_inDynamic"):
    #     if self._userModel()._inDynamic() or self._inUpdateWeight():
    #       newName = generateNameT(name, self._userModel().currentTimeStep())

    #   path = os.path.join(directoryPrefix, newName)
    #   # import pcraster
    #   # pcraster.report(variable, path)

    # def _readmapNew(self, name,
    #   style=1):
    #   """

    #   .. todo::

    #     `style` argument is not used.
    #   """
    #   directoryPrefix = ""
    #   nameSuffix = ".map"
    #   newName = ""

    #   if hasattr(self._userModel(), "_inStochastic"):
    #     if self._userModel()._inStochastic():
    #       if self._userModel()._inPremc() or self._userModel()._inPostmc():
    #         newName = name + nameSuffix
    #       else:
    #         directoryPrefix = str(self._userModel().currentSampleNumber())

    #   if hasattr(self._userModel(), "_inInitial"):
    #     if self._userModel()._inInitial():
    #       newName = name + nameSuffix

    #   if self._inResume():
    #     timestep = self._userModel().firstTimeStep()
    #     newName = generateNameT(name, timestep - 1)

    #   if hasattr(self._userModel(), "_inDynamic"):
    #     if self._userModel()._inDynamic() or self._inUpdateWeight():
    #       timestep = self._userModel().currentTimeStep()
    #       newName = generateNameT(name, timestep)

    #   path = os.path.join(directoryPrefix, newName)
    #   assert path is not ""
    #   # import pcraster
    #   # return pcraster.readmap(path)

    # def _assertAndThrow(self,
    #   expression,
    #   message):
    #   assert expression, message

class MonteCarloBase(object):
    def __init__(self):
        if self.__class__ is MonteCarloBase:
            raise NotImplementedError

        self._d_firstSampleNumber = 0
        self._d_lastSampleNumber = 0
        self._d_currentSampleNumber = 1
        self._d_inSample = False
        self._d_inStochastic = True
        self._d_inPremc = False
        self._d_inPostmc = False

    def premcloop(self):
        msg = "Class needs to implement 'premcloop' method"
        raise NotImplementedError(msg)

    def postmcloop(self):
        msg = "Class needs to implement 'postmcloop' method"
        raise NotImplementedError(msg)

    def nrSamples(self):
        """Return the number of samples."""
        assert self._d_firstSampleNumber
        return self._d_lastSampleNumber - self._d_firstSampleNumber + 1

    def currentSampleNumber(self):
        """Returns the current sample number."""
        assert self._d_currentSampleNumber
        return self._d_currentSampleNumber

    def sampleNumbers(self):
        """Returns a list of sample numbers configured."""
        assert self._d_firstSampleNumber
        return range(self._d_firstSampleNumber, self._d_lastSampleNumber + 1)

  def _inStochastic(self):
    if not hasattr(self, "_d_inStochastic"):
      return False
    return self._d_inStochastic

  def _inPremc(self):
    return self._d_inPremc

  def _inPostmc(self):
    return self._d_inPostmc

  def _lastSampleNumber(self):
    return self._d_lastSampleNumber

  def _firstSampleNumber(self):
    return self._d_firstSampleNumber

  def _setCurrentSample(self, nr):
    """
    Set the current sample number to nr.
    """
    assert nr >= self._firstSampleNumber()
    assert nr <= self._lastSampleNumber()
    self._d_currentSampleNumber = nr

  def _inSample(self):
    """
    Return whether a sample is currently executing.
    """
    #if hasattr(self._userModel(), "_d_inSample"):
    return self._d_inSample
    #else:
    #  return False

#  # -*- coding: utf-8 -*-
# import os
# import shutil
# import multiprocessing
# import sys
# import time
# from . import dynamicFramework
# from . import forkscript
# from . import frameworkBase
# from . import staticFramework

class MonteCarloFramework(frameworkBase.FrameworkBase, forkscript.ForkScript):
  """
  Framework class for the Monte Carlo method.

  `userModel`
    Instance that models the :ref:`Monte Carlo Model Concept <monteCarloModelConcept>`.

  `nrSamples`
    Number of realisations to run.

  `remove_dirs`
    Flag whether sample directories should be removed
  """

  def __init__(self,
    userModel,
    nrSamples=0,
    remove_dirs=True):
    frameworkBase.FrameworkBase.__init__(self)
    forkscript.ForkScript.__init__(self)

    self._d_model = userModel
    self._testRequirements()

    # adding framework specific attributes and methods

    # keep this fttb
    self._addMethodToClass(self._readmapNew)
    self._addMethodToClass(self._reportNew)

    assert isinstance(nrSamples, int)
    try:
      self._setSampleNumbers(1, nrSamples)
    except Exception as msg:
      sys.stderr.write('Error: %s\n' % str(msg))
      sys.exit(1)

    # Consecutive model runs by default
    self._d_forkSamples = False

    self._initialiseSampleDirectories(remove_dirs)

  def setForkSamples(self,
    fork,
    nrCPUs=1):
    """
    Set the forking of samples on or off.

    `fork`
      True or False.

    `nrCPUs`
      Number of CPU's to use. If not provided, this is autodectected.

    When forking is on EVERY sample will be forked to its own
    process. This is mainly useful on a cluster with automatic
    process migration or on SMP machines.

    .. warning::

       setrandomseed does not work when forking is enabled

    .. note::

       Support for forking samples is not available on Windows

    .. todo::

       Implement a portable solution for forking.

    .. todo::

       Default nrCPUs should be <= 0. The user may want to use only 1 CPU (yes,
       silly).
    """
    if (sys.platform == "win32"):
      self.showWarning("Forking not available on Windows platforms")
      self._d_nrProcessors = 1
      self._d_forkSamples = False
    else:
      # default detect nrCPUs
      if(nrCPUs == 1):
        self._d_nrProcessors = multiprocessing.cpu_count()
      else:
        self._d_nrProcessors = nrCPUs
      self._d_forkSamples = fork

  def _forkSamples(self):
    """
    Return whether samples should be forked.
    """
    return self._d_forkSamples

  def _userModel(self):
    """
    Return the class provided by the user.
    """
    return self._d_model._userModel()

  def _userFramework(self):
    """
    Return the framework provided by the user.
    """
    return self._d_model

  def _testRequirements(self):
    """
    Test whether the user model models the
    :ref:`Monte Carlo Model Concept <monteCarloModelConcept>`.

    .. todo::

       The implementation uses the staticFramework and dynamicFramework modules,
       which makes this class dependent on these specific framework classes.
       This is not necessary. The implementation should just perform the
       concept check and be done with it. Only use the interface defined by
       the concept, and nothing more. Don't assume anything about the
       implementation.
    """
    if (not isinstance(self._userFramework(), staticFramework.StaticFramework)) and (not isinstance(self._userFramework(), dynamicFramework.DynamicFramework)):
      msg = "Cannot run MonteCarlo framework: User model must be type of \
StaticFramework or DynamicFramework"
      raise frameworkBase.FrameworkError(msg)

    if not hasattr(self._userModel(), "premcloop"):
      self.showWarning("No premcloop section defined.")

    if not hasattr(self._userModel(), "postmcloop"):
      self.showWarning("No postmcloop section defined.")

  ## \brief Re-implemented from ShellScript.
  #
  # Runs the user model in the Monte Carlo mode.
  def run(self, premc=True, postmc=True):
    self._atStartOfScript()
    self._check()

    if premc:
      self._runPremcloop()

    # Samples.
    sample = self._userModel()._firstSampleNumber()
    while sample <= self._userModel()._lastSampleNumber():

      if self._forkSamples():
        assert self.isParentProcess()
        while self._systemIsOccupied(self.nrChildProcesses()):
          if self._debug():
            self.showMessage("System is occupied")
          # Wait for some seconds.
          # But first handle finished children.
          for childProcess in self.handleFinishedChildProcesses():
            if not self._quiet():
              if self._debug():
                self.showMessage(childProcess.message())
          time.sleep(0.1)

        # Throw in another sample.
        if self._debug():
          self.showMessage("Starting sample %d" % (sample))
        self.fork(str(sample))

        if self.isParentProcess():
          # Stop parent process to enable child process to start and
          # the load to adjust.
          # But first handle finished children.
          for childProcess in self.handleFinishedChildProcesses():
            if not self._quiet():
              if self._debug():
                self.showMessage(childProcess.message())
          time.sleep(0.01)

      # Child processes start here.
      if not self._forkSamples() or \
         (self._forkSamples() and self.isChildProcess()):
        self._incrementIndentLevel()
        self._atStartOfSample(sample)
        self._userModel()._setCurrentSample(sample)
        # Execute model
        self._d_model.run()

        self._sampleFinished()
        self._decrementIndentLevel()

      # Child processes stop here.
      if self._forkSamples():
        if self.isChildProcess():
          os._exit(0)

      assert not self._forkSamples() or self.isParentProcess()
      sample += 1

    if self._forkSamples() and  self.isParentProcess():
      for childProcess in self.waitForChildProcessesToFinish():
        if not self._quiet():
          if self._debug():
            self.showMessage(childProcess.message())
      assert not self.childProcesses()

    if postmc:
      self._runPostmcloop()

    return 0



  ## \brief Creates the directories in which the sample data can be stored.
  #
  # \attention Already existing sample directories will be cleaned!
  def _initialiseSampleDirectories(self, remove_sample_dirs):
    sample = self._userModel()._firstSampleNumber()
    while sample <= self._userModel()._lastSampleNumber():
      dirname = '{}'.format(sample)

      if not os.path.exists(dirname):
        # Create sample directory.
        os.mkdir(dirname)
      elif remove_sample_dirs == True:
        if not os.path.isdir(dirname):
          # Remove existing file with name of sample directory.
          os.remove(dirname)
          os.mkdir(dirname)
        else:
          # remove existing and create emtpy directories
          shutil.rmtree(dirname)
          os.mkdir(dirname)

      assert os.path.exists(dirname) and os.path.isdir(dirname)
      sample += 1


  def _setSampleNumbers(self, firstSampleNumber, lastSampleNumber):
    assert firstSampleNumber > 0
    if lastSampleNumber <= 0:
      msg = "number of samples argument (%s) of MonteCarloFramework must be > 0" % (lastSampleNumber)
      raise AttributeError(msg)
    assert lastSampleNumber >= firstSampleNumber

    self._userModel()._d_firstSampleNumber = firstSampleNumber
    self._userModel()._d_lastSampleNumber = lastSampleNumber

  ## \brief Returns true if new processes can be forked
  def _systemIsOccupied(self, nrChilds):
    return nrChilds >= self._d_nrProcessors

  ## Checks the current configuration of the script.
  def _check(self):
    if self._userModel().nrSamples() == 0:
      self.showWarning("""
   Since the number of samples to execute is 0, only the premcloop
   and the postmcloop functions will be executed. Any script that does
   something meaningful needs to calculate at least one "sample".""")

  def setQuiet(self, quiet=True):
    """
    Disables the progress display of sample numbers.
    """
    self._d_quietProgressSampleNr = quiet


class DynamicFramework(FrameworkBase):
    """
    Framework class for dynamic models.

    `userModel`
      Instance that models the :ref:`Dynamic Model Concept <dynamicModelConcept>`.

    `lastTimeStep`
      Last timestep to run.

    `firstTimestep`
      Sets the starting timestep of the model (optional, default is 1).
    """

    def __init__(self,
                 userModel,
                 lastTimeStep=0,
                 firstTimestep=1):

        FrameworkBase.__init__(self)
        self._d_model = userModel
        self._testRequirements()
        # # fttb
        # self._addMethodToClass(self._readmapNew)
        # self._addMethodToClass(self._reportNew)
        try:
            self._userModel()._setNrTimeSteps(lastTimeStep)
            self._d_firstTimestep = firstTimestep
            self._userModel()._setFirstTimeStep(self._d_firstTimestep)
        except Exception as msg:
            sys.stderr.write('Error: %s\n' % str(msg))
            sys.exit(1)

    def _userModel(self):
        """
        Return the model instance provided by the user.
        """
        return self._d_model

    def run(self):
        """Run the dynamic user model.
        .. todo::
        This method depends on the filter frameworks concept. Shouldn't its run
        method call _runSuspend()?
        """
        self._atStartOfScript()
        if(hasattr(self._userModel(), "resume")):
            if self._userModel().firstTimeStep() == 1:
                self._runInitial()
            else:
                self._runResume()
        else:
            self._runInitial()

        self._runDynamic()
        # Only execute this section while running filter frameworks.
        if hasattr(self._userModel(), "suspend") and \
           hasattr(self._userModel(), "filterPeriod"):
            self._runSuspend()
        return 0

    def _testRequirements(self):
        pass
        # """
        # Test whether the user model models the
        # :ref:`Dynamic Model Concept <dynamicModelConcept>`.
        # """
        # if hasattr(self._userModel(), "_userModel"):
        #     msg = "The _userModel method is deprecated and obsolete"
        #     self.showWarning(msg)

        # if( not hasattr(self._userModel(), "dynamic") and not hasattr(self._userModel(), "run")):
        #     msg = "Cannot run dynamic framework: Implement dynamic method"
        #     raise FrameworkError(msg)

        # if not hasattr(self._userModel(), "initial"):
        #     if self._debug():
        #         self.showWarning("No initial section defined.")

    def setQuiet(self, quiet=True):
        self._d_quietProgressDots = quiet

class HmDynamicFramework(DynamicFramework):
    pass

class HmMonteCarloFramework(MonteCarloFramework):
    pass

class HmEnsKalmanFilterFramework(EnsKalmanFilterFramework):
    pass

