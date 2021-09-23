#!/usr/bin/env python3

import os
import sys
import shutil
import time
import multiprocessing
# from . import dynamicFramework
# from . import staticFramework
from .forkscript import ForkScript
from .frameworkbase import FrameworkBase

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
        """Set the current sample number to nr."""
        assert nr >= self._firstSampleNumber()
        assert nr <= self._lastSampleNumber()
        self._d_currentSampleNumber = nr

    def _inSample(self):
        """Return whether a sample is currently executing."""
        return self._d_inSample


class MonteCarloFramework(FrameworkBase, ForkScript):
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
                 frameworkbase.FrameworkBase.__init__(self)
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

    def setForkSamples(self, fork, nrCPUs=1):
        """Set the forking of samples on or off.
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
           Default nrCPUs should be <= 0. The user may want to use only 1 CPU
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
        """Return whether samples should be forked."""
        return self._d_forkSamples

    def _userModel(self):
        """Return the class provided by the user."""
        return self._d_model._userModel()

    def _userFramework(self):
        """Return the framework provided by the user."""
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
          self.showWarning(
              """Since the number of samples to execute is 0, only the premcloop
              and the postmcloop functions will be executed. Any script that does
              something meaningful needs to calculate at least one "sample".
              """
          )

    def setQuiet(self, quiet=True):
        """Disables the progress display of sample numbers."""
        self._d_quietProgressSampleNr = quiet
