#!/usr/bin/env python
# -*- coding: utf-8 -*-

from configparser import ConfigParser, ExtendedInterpolation
import os
import time
import datetime
import shutil
import glob

# from . import file_handling
# from . import disclaimer
# from .Messages import ModelError, ModelFileError, ModelWarning

import logging
logger = logging.getLogger(__name__)

class Configuration(object):

    def __init__(
            self,
            config_filename,
            debug_mode = False# ,
            # no_modification = True,
            # system_arguments = None,
            # relative_ini_weather_paths = False
    ):
        object.__init__(self)
        if config_filename is None:
            raise ValueError('No configuration file specified')

        self._timestamp = datetime.datetime.now()
        self.iniFileName = os.path.abspath(iniFileName)
        self.debug_mode = debug_mode
        self.parse_configuration_file(self.iniFileName)
        if no_modification: self.set_configuration(system_arguments)

    def set_configuration(self, system_arguments = None):
        self.repair_ini_key_names()
        self.set_clone_map()
        self.set_land_mask()
        self.create_output_directories()
        self.create_coupling_directories()        
        self.repair_logging_key_names()
        self.initialize_logging("Default", system_arguments)
        self.backup_configuration()

    def initialize_logging(self, log_file_location = "Default", system_arguments = None):
        """Initialize logging. Prints to both the console 
        and a log file, at configurable levels.
        """
        logging.getLogger().setLevel(logging.DEBUG)
        formatter = logging.Formatter('%(asctime)s %(name)s %(levelname)s %(message)s')
        log_level_console = self.LOGGING['log_level_console']        
        log_level_file = self.LOGGING['log_level_file']        

        # log level for debug mode:
        if self.debug_mode == True: 
            log_level_console = "DEBUG"
            log_level_file    = "DEBUG"

        console_level = getattr(logging, log_level_console.upper(), logging.INFO)
        if not isinstance(console_level, int):
            raise ValueError('Invalid log level: %s', log_level_console)
        
        # create handler, add to root logger
        console_handler = logging.StreamHandler()
        console_handler.setFormatter(formatter)
        console_handler.setLevel(console_level)
        logging.getLogger().addHandler(console_handler)

        # log file name (and location)
        if log_file_location != "Default":  self.logFileDir = log_file_location
        log_filename = os.path.join(
            self.logFileDir,
            os.path.basename(self.iniFileName) + '_'
            + str(self._timestamp.isoformat()).replace(":",".")
            + '.log'
        )

        file_level = getattr(logging, log_level_file.upper(), logging.DEBUG)
        if not isinstance(console_level, int):
            raise ValueError('Invalid log level: %s', log_level_file)

        # create handler, add to root logger
        file_handler = logging.FileHandler(log_filename)
        file_handler.setFormatter(formatter)
        file_handler.setLevel(file_level)
        logging.getLogger().addHandler(file_handler)
        
        dbg_filename = os.path.join(
            self.logFileDir,
            os.path.basename(self.iniFileName) + '_'
            +  str(self._timestamp.isoformat()).replace(":",".")
            + '.dbg'
        )
        
        # create handler, add to root logger
        debug_handler = logging.FileHandler(dbg_filename)
        debug_handler.setFormatter(formatter)
        debug_handler.setLevel(logging.DEBUG)
        logging.getLogger().addHandler(debug_handler)

        disclaimer.print_disclaimer(with_logger = True)        
        logger.info('Model run started at %s', self._timestamp)
        logger.info('Logging output to %s', log_filename)
        logger.info('Debugging output to %s', dbg_filename)
        
        if system_arguments != None:
            logger.info('The system arguments given to execute this run: %s', system_arguments)        

    def backup_configuration(self):
        """Function to copy ini file to log directory"""
        backup_location = os.path.join(
            self.logFileDir,
            os.path.basename(self.iniFileName) + '_'
            + str(self._timestamp.isoformat()).replace(":",".")
            + '.ini'
        )
        shutil.copy(self.iniFileName, backup_location)
        
    def parse_configuration_file(self, modelFileName):
        config = ConfigParser(interpolation=ExtendedInterpolation())
        config.optionxform = str
        config.read(modelFileName)
        self.allSections = config.sections()
        for section in self.allSections:
            vars(self)[section] = {}
            options = config.options(section)
            for option in options:
                val = config.get(section, option)
                self.__getattribute__(section)[option] = val

    def set_clone_map(self):
        try:
            self.cloneMap = str(self.MODEL_GRID['cloneMap'])
        except:
            self.cloneMap = None

    def set_land_mask(self):
        try:
            self.landmask = str(self.MODEL_GRID['landmask'])
        except:
            self.landmask = None
            
    def create_output_directories(self):
        cleanOutputDir = False
        if cleanOutputDir:
            try:
                shutil.rmtree(self.FILE_PATHS['PathOut'])
            except:
                pass

        try:
            os.makedirs(self.FILE_PATHS['PathOut'])
        except:
            pass

        self.tmpDir = os.path.join(self.FILE_PATHS['PathOut'], 'tmp')
        if os.path.exists(self.tmpDir):
            shutil.rmtree(self.tmpDir)
        os.makedirs(self.tmpDir)        

        self.outNCDir = os.path.join(self.FILE_PATHS['PathOut'], 'netcdf')
        if os.path.exists(self.outNCDir):
            shutil.rmtree(self.outNCDir)
        os.makedirs(self.outNCDir)

        # self.scriptDir = os.path.join(self.FILE_PATHS['PathOut'], 'scripts')
        # if os.path.exists(self.scriptDir):
        #     shutil.rmtree(self.scriptDir)
        # os.makedirs(self.scriptDir)

        # # copy scripts from the current directory to the backup 
        # path_of_this_module = os.path.abspath(os.path.dirname(__file__))
        # self.starting_directory = path_of_this_module
        # all_files = glob.glob(os.path.join(path_of_this_module, '*.py'))
        # for filename in all_files:
        #     shutil.copy(filename, self.scriptDir)
        
        self.logFileDir = os.path.join(self.FILE_PATHS['PathOut'], 'log')        
        cleanLogDir = True
        if os.path.exists(self.logFileDir) and cleanLogDir:
            shutil.rmtree(self.logFileDir)
        os.makedirs(self.logFileDir)

        self.endStateDir = os.path.join(self.FILE_PATHS['PathOut'], 'states')
        if os.path.exists(self.endStateDir):
            shutil.rmtree(self.endStateDir)
        os.makedirs(self.endStateDir)

    def create_coupling_directories(self):
        pass

    def repair_logging_key_names(self):
        if 'LOGGING' not in self.allSections:
            self.LOGGING = {}
        self.LOGGING['log_level_console'] = "INFO"
        self.LOGGING['log_level_file'] = "INFO"
        
    def repair_ini_key_names(self):
        """This function is used to change or modify key names of
        options, to check the validity of options and to infill 
        missing keys"""
        pass
