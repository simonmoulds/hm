#!/usr/bin/env python3

# import atexit
# import os
# import re
import sys
# import weakref

# TODO: Remove this error class
class FrameworkError(Exception):
  def __init__(self,
    msg):
    self._msg = msg

  def __str__(self):
    return self._msg


# TODO look up the atexit module
# @atexit.register
# def _atExit():
#   print


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
