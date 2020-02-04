#!/usr/bin/env python
# -*- coding: utf-8 -*-

from configparser import ConfigParser, ExtendedInterpolation
import os
import time
# import datetime
import pandas as pd
import shutil
import glob

# from . import file_handling
# from . import disclaimer
# from .Messages import ModelError, ModelFileError, ModelWarning

import logging
logger = logging.getLogger(__name__)

required_config_sections = []

class Configuration(object):
    def __init__(
            self,
            config_filename,
            debug_mode = False
    ):
        if config_filename is None:
            raise ValueError(
                'No configuration file specified'
            )
        
        config_filename = os.path.abspath(config_filename)
        if not os.path.exists(config_filename):
            raise ValueError(
                'Specified configuration file '
                + config_filename + ' does not exist'
            )
        self.config_filename = config_filename
        self._timestamp = pd.Timestamp.now()
        self._debug_mode = debug_mode
        self.parse_config_file(self.config_filename)
        self.set_config(system_arguments)

    def parse_config_file(self):
        """Parse the configuration file."""
        config = ConfigParser(interpolation=ExtendedInterpolation())
        config.optionxform = str
        config.read(self.config_filename)
        self.config_sections = config.sections()
        for section in self.config_sections:
            vars(self)[section] = {}
            options = config.options(section)
            for option in options:
                val = config.get(section, option)
                self.__getattribute__(section)[option] = val
        
    def set_config(self, system_arguments = None):
        self.check_required_options()
        self.assign_default_options()
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
        if self._debug_mode == True: 
            log_level_console = "DEBUG"
            log_level_file    = "DEBUG"

        console_level = getattr(logging, log_level_console.upper(), logging.INFO)
        if not isinstance(console_level, int):
            raise ValueError('Invalid log level: %s', log_level_console)
        
        console_handler = logging.StreamHandler()
        console_handler.setFormatter(formatter)
        console_handler.setLevel(console_level)
        logging.getLogger().addHandler(console_handler)

        if log_file_location != "Default":
            self.logFileDir = log_file_location
        log_filename = os.path.join(
            self.logFileDir,
            os.path.basename(self.config_filename) + '_'
            + str(self._timestamp.isoformat()).replace(":",".")
            + '.log'
        )
        file_level = getattr(logging, log_level_file.upper(), logging.DEBUG)
        if not isinstance(console_level, int):
            raise ValueError('Invalid log level: %s', log_level_file)

        file_handler = logging.FileHandler(log_filename)
        file_handler.setFormatter(formatter)
        file_handler.setLevel(file_level)
        logging.getLogger().addHandler(file_handler)
        
        dbg_filename = os.path.join(
            self.logFileDir,
            os.path.basename(self.config_filename) + '_'
            +  str(self._timestamp.isoformat()).replace(":",".")
            + '.dbg'
        )
        
        debug_handler = logging.FileHandler(dbg_filename)
        debug_handler.setFormatter(formatter)
        debug_handler.setLevel(logging.DEBUG)
        logging.getLogger().addHandler(debug_handler)

        disclaimer.print_disclaimer(with_logger = True)        
        logger.info('Model run started at %s', self._timestamp)
        logger.info('Logging output to %s', log_filename)
        logger.info('Debugging output to %s', dbg_filename)
        
        if system_arguments != None:
            logger.info(
                'The system arguments given to execute this run: %s',
                system_arguments
            )        

    def backup_configuration(self):
        """Copy config file to log directory."""
        backup_location = os.path.join(
            self.logFileDir,
            os.path.basename(self.iniFileName) + '_'
            + str(self._timestamp.isoformat()).replace(":",".")
            + '.ini'
        )
        shutil.copy(self.iniFileName, backup_location)
        
    # def set_clone_map(self):
    #     try:
    #         self.cloneMap = str(self.MODEL_GRID['cloneMap'])
    #     except:
    #         self.cloneMap = None

    # def set_land_mask(self):
    #     try:
    #         self.landmask = str(self.MODEL_GRID['landmask'])
    #     except:
    #         self.landmask = None            
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

    def assign_default_options(self):
        self._assign_default_logging_options()
        self._assign_default_file_path_options()
                    
    def assign_default_logging_options(self):
        if 'LOGGING' not in self.config_sections:
            self.LOGGING = {}
        if 'log_level_console' not in self.LOGGING:
            self.LOGGING['log_level_console'] = 'INFO',
        if 'log_level_file' not in self.LOGGING:
            self.LOGGING['log_level_file'] = 'INFO',
            
    def assign_default_file_path_options(self):
        if 'FILE_PATHS' not in self.config_sections:
            self.FILE_PATHS = {}
        if 'in_path' not in self.FILE_PATHS:
            self.FILE_PATHS['in_path'] = '.'
        if 'out_path' not in self.FILE_PATHS:
            self.FILE_PATHS['in_path'] = '.'

    def assign_default_reporting_options(self):
        if 'REPORTING' not in self.config_sections:
            self.REPORTING = {}
        if 'report' not in self.REPORTING:
            self.REPORTING['report'] = False
            
    def check_required_options(self):
        self._default_config_check('MODEL_GRID', ['mask'])
        
    def _default_config_check(self, section, options):
        if section not in self.config_sections:
            raise KeyError(
                self.generate_missing_section_message(section)
            )
        else:
            for option in options:
                if option not in vars(self)[section]:
                    raise KeyError(
                        self.generate_missing_option_message(section, option)
                    )
                
    def generate_missing_section_message(self, section):
        return 'Configuration file ' + self.config_filename \
            + ' does not contain section ' + section \
            + ', which must be supplied'

    def generate_missing_option_message(self, section, option):
        return 'Section ' + section + ' in configuration file ' \
            + self.config_filename + ' does not contain option ' \
            + option + ', which must be supplied'


