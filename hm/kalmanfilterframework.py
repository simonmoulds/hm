#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
# import random
import sys
import shutil
import numpy
import pickle
from numpy import linalg
# from . import dynamicframework
from . import frameworkbase
from . import montecarloframework
# from .frameworkbase import generateNameT#, generateNameS, generateNameST

class EnKfBase(object):
    def __init__(self):
        if self.__class__ is EnKfBase:
            raise NotImplementedError

    def initial(self):
        msg = "Class needs to implement 'initial' method"
        raise NotImplementedError(msg)

    def setDebug(self):
        msg = "Class needs to implement 'setDebug' method"
        raise NotImplementedError(msg)

    def setState(self):
        msg = "Class needs to implement 'setState' method"
        raise NotImplementedError(msg)

    def setObservations(self):
        msg = "Class needs to implement 'setObservations' method"
        raise NotImplementedError(msg)


class EnKfModel(EnKfBase):
    def __init__(self):
        EnKfBase.__init__(self)
    ## \brief Storing map data to disk.
    #
    # \param variable object containing the PCRaster map data
    # \param name name used as filename. Use filename without extension. File extension for deterministic
    # static models is ".map" and will be appended automatically.
    def report(self, variable, name, style=1):
        self._reportNew(variable, name)

    ## \brief Reads map data from disk.
    #
    # \param name name used as filename. Use filename without extension. File extension for deterministic
    # static models is ".map" and will be appended automatically.
    def readmap(self, name, style=1):
        return self._readmapNew(name)


class EnsKalmanFilterFramework(frameworkbase.FrameworkBase):
    def __init__(self, userModel):
        frameworkbase.FrameworkBase.__init__(self)
        self._d_model = userModel
        self._testRequirements()
        self._d_totalTimesteps = self._userModel().nrTimeSteps()
        self._d_trackCloned = {}
        # adding framework specific attributes and methods
        self._addAttributeToClass("_d_filterPeriod", 0)
        self._addAttributeToClass("_d_inFilterPeriod", False)
        self._addAttributeToClass("_d_filterTimesteps", [])
        self._addAttributeToClass("_d_inResume", False)
        self._addAttributeToClass("_d_inUpdateWeight", False)
        self._resetSampleWeights()

        self._addMethodToClass(self.getStateVector)
        self._addMethodToClass(self._runPremcloop)
        self._addMethodToClass(self._runPostmcloop)

        # self._addMethodToClass(self.readmap)
        # self._addMethodToClass(self.readDeterministic)
        self._addMethodToClass(self.setMeasurementOperator)
        self._addMethodToClass(self.setObservedMatrices)

        # \todo !!!test if filter timesteps are in interval of model timesteps...
        self.sizeStateVector = 0
        self._initialiseObservedDir()

    def setObservedMatrices(self, observations, covariance):
        assert type(observations) == numpy.ndarray
        assert type(covariance) == numpy.ndarray
        filtermoment = self._userModel().currentTimeStep()

        fileName = os.path.join("observedState",'obs%s.tmp' % (filtermoment))
        file = open(fileName, 'wb')
        pickle.dump(observations, file)
        file.close()
        fileName = os.path.join("observedState",'cov%s.tmp' % (filtermoment))
        file = open(fileName, 'wb')
        pickle.dump(covariance, file)
        file.close()

    def setMeasurementOperator(self, matrix):
        # If this is not used the identify matrix will be used
        assert type(matrix) == numpy.ndarray
        filtermoment = self._userModel().currentTimeStep()
        fileName = os.path.join("observedState",'h%s.tmp' % (filtermoment))
        file = open(fileName, 'wb')
        pickle.dump(matrix, file)
        file.close()

    def _testRequirements(self):
        if not isinstance(self._d_model, montecarloframework.MonteCarloFramework):
            self.showError("Model must be instance of MonteCarloFramework.")
            sys.exit()

        if not hasattr(self._d_model, 'run'):
            self.showError("No 'run' section defined.")
            sys.exit()

        if not hasattr(self._userModel(), 'setState'):
            self.showError("No 'setState' function defined.")
            sys.exit()

        if not hasattr(self._userModel(), 'resume'):
            msg = "Cannot run particle filter framework: Implement 'resume' method"
            raise frameworkBase.FrameworkError(msg)

    def _particleWeights(self):
        return self._userModel()._d_particleWeights

    def _userModel(self):
        return self._d_model._userModel()

    def _initialiseObservedDir(self):
        varName = "observedState"
        if not os.path.isdir(varName):
            os.mkdir(varName)
        else:
            shutil.rmtree(varName)
            os.mkdir(varName)

    def _initialiseStateDir(self):
        # Create subdirectories for state variables
        varName = "stateVector"
        if not os.path.isdir(varName):
            # Create sample directory.
            os.mkdir(varName)
        else :
            #if not os.path.isdir(varName):
            #  # Remove existing file with name of sample directory.
            shutil.rmtree(varName)
            os.mkdir(varName)

    def _initialiseSampleDirectories(self):
        sample = self._userModel()._firstSampleNumber()
        while sample <= self._userModel()._lastSampleNumber():
            cwd = os.getcwd()
            dirname = "%d" % (sample)
            varName = "stateVar"
            os.chdir(dirname)
            if not os.path.isdir(varName):
                os.mkdir(varName)
            else:
                os.remove(varName)
                os.mkdir(varName)
            os.chdir(cwd)
            assert os.path.exists(os.path.join(dirname,"stateVar")) and os.path.isdir(os.path.join(dirname,"stateVar"))
            sample += 1

    def setFilterTimesteps(self, filterTimesteps):
        # Set the filter moments
        assert type(filterTimesteps) == list or type(filterTimesteps) == numpy.ndarray
        for filtertimestep in filterTimesteps:
            assert filtertimestep < self._userModel().nrTimeSteps()
        self._userModel()._d_filterTimesteps = filterTimesteps

    def filterTimesteps(self):
        # Return a list of filter moments
        return self._userModel()._d_filterTimesteps

    def run(self):
        # Run the user model in filter mode
        if hasattr(self._userModel(), 'run'):
            self._userModel().run()
        else:
            self._atStartOfScript()
            self._initialiseStateDir()
            self._initialiseSampleDirectories()
            lastPeriod = len(self._userModel()._d_filterTimesteps)
            if lastPeriod == 0:
                self.showError("No filter timesteps specified")
                sys.exit()

            # set the proposal/initial weight distribution by user
            if hasattr(self._userModel(), 'setInitialParticleWeights'):
                self._userModel()._d_particleWeights = self._userModel().setInitialParticleWeights()
                # check initial weights
                assert type(self._particleWeights()) == list
                assert len(self._particleWeights()) == self._userModel().nrSamples()
                for i in range(0, len(self._particleWeights())):
                    assert type(self._particleWeights()[i]) == float

            self._userModel()._runPremcloop()

            # looping over the filter periods
            for currentPeriod in range(0, len(self._userModel()._d_filterTimesteps) + 1):
                # \todo replace with a better solution...
                sumW = sum(self._particleWeights())
                assert abs(sumW - 1.0) < 0.00001
                self._runMonteCarlo(currentPeriod, lastPeriod)

                if not currentPeriod == lastPeriod:
                    # retrieve the state vectors for each sample
                    for sample in range(1, self._userModel().nrSamples() + 1):
                        self._userModel()._setCurrentSample(sample)
                        self._userModel()._d_inUpdateWeight = True
                        stateVector = self._userModel().setState()
                        self._userModel()._d_inUpdateWeight = False
                        assert type(stateVector) == numpy.ndarray
                        fileName = os.path.join("stateVector",'ensMember%s.tmp' %(sample))
                        file = open(fileName,'wb')
                        pickle.dump(stateVector, file)
                        file.close()

                    # for current update moment
                    self._getObservedValues()
                    self._kalmanFilter()
                    currentPeriod += 1
                    self._userModel()._d_filterPeriod += 1

        self._userModel()._setFirstTimeStep(1)
        self._userModel()._runPostmcloop()
        return 0

    def _getObservedValues(self):
      self._userModel().setObservations()

    def _kalmanFilter(self):
        """Implementation of Ensemble Kalman Filter.

        This follows equations 44-52 from Geir Evensen's paper
        'The Ensemble Kalman Filter: theoretical formulation
        and practical implementation'

        n size of state vector (sizeStateVector)
        m nr of observations (sizeObservedVector)
        N nr of ensemble members

        A matrix with model states
        H matrix 'measurement operator'
        D matrix with observations
        """

        # ##############################
        # Load current model state
        # ##############################
        fileName = os.path.join("stateVector",'ensMember%s.tmp' %(str(1)))
        file = open(fileName,'rb')
        vec = pickle.load(file)
        sizeStateVector = len(vec)
        file.close()

        # ##############################
        # Load observations
        # ##############################

        # length of the observed vector \todo do we know that?
        fileName = os.path.join("observedState","obs%s.tmp" %(self._userModel()._d_filterTimesteps[self._userModel()._d_filterPeriod]))
        file = open(fileName,'rb')
        vec = pickle.load(file)
        sizeObservedVector = len(vec)
        file.close()

        # number of ensemble members
        nrEnsembleMembers =  self._userModel().nrSamples()

        # create A
        A = numpy.zeros((sizeStateVector, nrEnsembleMembers), dtype=float)

        # ##############################
        # Fill A with current model state
        # ##############################

        # NB - we seem to be loading stateVector twice?
        # \todo is there a better way to construct a matrix from vecors?
        for sample in range(1, self._userModel().nrSamples() + 1):
            fileName = os.path.join("stateVector",'ensMember%s.tmp' %(sample))
            file = open(fileName,'rb')
            vec = pickle.load(file)
            file.close()
            # A[:,sample-1] = vec
            for i in range(0, sizeStateVector):
                A[i,sample-1] = vec[i]

        # ##############################
        # Compute measurement operator
        # ##############################

        # H is the measurement operator which relates the true model
        # state to the observations, allowing for measurement errors

        fileName = os.path.join("observedState","h%s.tmp" %(self._userModel()._d_filterTimesteps[self._userModel()._d_filterPeriod]))
        if os.path.exists(fileName):
            # First, we attempt to obtain user-specified H
            file = open(fileName,'rb')
            H = pickle.load(file)
            file.close()
        else:
            # Otherwise use the identity matrix - this will convert
            # the state to measurement 1:1
            H = numpy.eye(sizeObservedVector, sizeStateVector, dtype=float)
        assert H.shape == (sizeObservedVector, sizeStateVector), "Shape of provided matrix H %s does not match (%s, %s)" %(H.shape, sizeObservedVector, sizeStateVector)

        # ##############################
        # Obtain measurements (D)
        # ##############################

        fileName = os.path.join("observedState","obs%s.tmp" %(self._userModel()._d_filterTimesteps[self._userModel()._d_filterPeriod]))
        file = open(fileName, 'rb')
        D = pickle.load(file)
        file.close()
        assert D.shape == (sizeObservedVector, nrEnsembleMembers), "Shape of provided matrix D %s does not match (%s, %s)" %(D.shape, sizeObservedVector, nrEnsembleMembers)

        # ##############################
        # Obtain the measurement error covariance matrix (Re)
        # ##############################

        # Equation 51

        fileName = os.path.join("observedState","cov%s.tmp" %(self._userModel()._d_filterTimesteps[self._userModel()._d_filterPeriod]))
        file = open(fileName, 'rb')
        Re = pickle.load(file)
        file.close()
        assert Re.shape == (sizeObservedVector, sizeObservedVector), "Shape of provided matrix Re %s does not match (%s, %s)" %(Re.shape, sizeObservedVector, sizeObservedVector)

        # ##############################
        # Calculate ensemble covariance matrix (Pe)
        # ##############################

        # Abar is simply the mean of each row, where each row represents
        # a vector of observed values across the ensemble, and each
        # column represents a vector of observations from the same
        # ensemble member (e.g. different values in space)

        # Equation 45
        Abar = numpy.dot(A, numpy.array( [[1.0/nrEnsembleMembers] * nrEnsembleMembers ] * nrEnsembleMembers, dtype=float))
        # Equation 46
        Ad = A - Abar
        # Equation 47
        Pe =  1.0/(nrEnsembleMembers - 1) * numpy.dot(Ad,numpy.transpose(Ad))

        # ##############################
        # Update A matrix
        # ##############################

        # the following variables reference terms in Equation 52
        DmAH = D - numpy.dot(H,A)                # D - HA
        PeHt = numpy.dot(Pe,numpy.transpose(H))  # Pe * H^T [^T denotes transpose]
        HPeHt = numpy.dot(H, PeHt)               # H * Pe * H^T
        HPeHtpRe = HPeHt + Re                    # (H * Pe * H^T + Re)
        INV = linalg.pinv(HPeHtpRe)              # ...
        INVDmAH = numpy.dot(INV, DmAH)           # ...
        A = A + numpy.dot(PeHt, INVDmAH)         # Equation 52
        for sample in range(1, self._userModel().nrSamples() + 1):
            fileName = os.path.join("stateVector",'a%s.tmp' %(sample))
            file = open(fileName,'wb')
            index = sample - 1
            vec = A[:,index]
            pickle.dump(vec, file)
            file.close()

    def getStateVector(self, sampleNumber):
        # Returns the updated variables
        fileName = os.path.join("stateVector",'a%s.tmp' %(sampleNumber))
        file = open(fileName,'rb')
        vec = pickle.load(file)
        file.close()
        return vec

    def _normaliseWeights(self, weights):
        assert weights
        sumWeights = sum(weights)
        norm = [0.0] * len(weights)
        for i in range(0, len(weights)):
            norm[i] =  weights[i] / sumWeights
        return norm

    def _resetSampleWeights(self):
        assert self._userModel().nrSamples() > 0
        self._userModel()._d_particleWeights = [1.0 / self._userModel().nrSamples()] * self._userModel().nrSamples()

    def _cumulativeWeights(self, weights):
        cumulative = [0.0] * self._userModel().nrSamples()
        value = 0.0
        for i in range(len(weights)):
            value += weights[i]
            cumulative[i] = value
        return cumulative

    def _startEndOfPeriod(self, currentPeriod, lastPeriod):
        # Determine start end end timestep of current period
        if currentPeriod == 0:
            startTimestep = 1
            endTimestep = self._userModel()._d_filterTimesteps[currentPeriod]
        elif currentPeriod == lastPeriod:
            startTimestep = self._userModel()._d_filterTimesteps[currentPeriod -1] + 1
            endTimestep = self._d_totalTimesteps
        else:
            startTimestep = self._userModel()._d_filterTimesteps[currentPeriod - 1] + 1
            endTimestep = self._userModel()._d_filterTimesteps[currentPeriod]
        assert startTimestep <= endTimestep
        return startTimestep, endTimestep

    def _executePrePostMc(self, currentPeriod, lastPeriod):
        if currentPeriod == 0:
            # execute premc
            premc = True
            postmc = False
        elif currentPeriod == lastPeriod:
            # execute postmc
            premc = False
            postmc = True
        else:
            # without pre/postmc
            premc = False
            postmc = False
        # TODO: assert something here?
        return premc, postmc

    def _runMonteCarlo(self, currentPeriod, lastPeriod):
        # Get user model and (re)set start and end time
        startTimestep, endTimestep = self._startEndOfPeriod(currentPeriod, lastPeriod)
        self._userModel()._setNrTimeSteps(endTimestep)
        self._userModel()._setFirstTimeStep(startTimestep)
        self._userModel()._setCurrentTimeStep(endTimestep)
        # Run the model in mc mode for current filter period
        self._incrementIndentLevel()
        self._atStartOfFilterPeriod(currentPeriod)
        self._d_model.run(False, False)
        self._atEndOfFilterPeriod()
        self._decrementIndentLevel()

    # # NOT USED
    # def readmap(self, name):
    #     # Read sample data from disk
    #     return self._readmapNew(name)

    # # NOT USED
    # def readDeterministic(self, name):
    #     # Read deterministic data from disk
    #     if self._userModel()._inPremc() or self._userModel()._inPostmc() or self._userModel()._inInitial():
    #         newName = name + ".map"
    #     else:
    #         newName = generateNameT(name, self._userModel().currentTimeStep())
