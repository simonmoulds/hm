#!/usr/bin/env python
# -*- coding: utf-8 -*-

# TODO: inherit from Basic Model Interface

class Model(object):
    def __init__(
            self,
            config,
            time,
            domain,
            is_1d,
            init=None,
            **kwargs
    ):
        """Model.
        
        Parameters
        ----------
        config: `hm.Configuration`
            Model configuration.
        time: `hm.ModelTime`
            Model time.
        domain: `hm.HmDomain`
            Model spatial domain.
        is_1d: bool
            Whether the model runs over a one-dimensional 
            representation of space.
        init: TODO
            Initial model state
        **kwargs: Any
            Not implemented.
        """
        self.config = config
        self.time = time
        self.domain = domain
        self.is_1d = is_1d
        self.init = init
        # EXPERIMENTAL:
        self.currentSampleNumber = 1
        self.filter_timesteps = []
