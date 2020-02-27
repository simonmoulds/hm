#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import shutil
import sys
import math
import gc
import numpy as np 
import rasterio
import xarray as xr

from .api import set_domain

class Model(object):
    def __init__(
            self,
            config,
            time,
            domain,
            init=None,
            **kwargs
    ):
        self.config = config
        self.time = time
        self.domain = domain
