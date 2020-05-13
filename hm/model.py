#!/usr/bin/env python
# -*- coding: utf-8 -*-


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
