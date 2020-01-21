#!/usr/bin/env python
# -*- coding: utf-8 -*-

import logging
logger = logging.getLogger(__name__)

def print_disclaimer(with_logger = False):

    disclaimer_message  =                                                                                                  "\n"
    disclaimer_message +=                                                                                                  "\n"
    disclaimer_message += " HydroModelBuilder                                                                          " + "\n"
    disclaimer_message +=                                                                                                  "\n"
    disclaimer_message += " Copyright (C) 2018, Simon Moulds, Wouter Buytaert,                                         " + "\n"
    disclaimer_message += " Dept. of Civil and Environmental Engineering, Imperial College London, United Kingdom      " + "\n"
    disclaimer_message +=                                                                                                  "\n"
    disclaimer_message += " This program comes with ABSOLUTELY NO WARRANTY                                             " + "\n"
    disclaimer_message += " This is free software, and you are welcome to redistribute it under certain conditions     " + "\n"
    disclaimer_message += " See the LICENSE file for more details                                                      " + "\n"
    disclaimer_message +=                                                                                                  "\n"
    
    if with_logger:
        logger.info(disclaimer_message)
    else:
        print(disclaimer_message)
