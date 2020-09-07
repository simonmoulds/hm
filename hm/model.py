#!/usr/bin/env python
# -*- coding: utf-8 -*-


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
        self.config = config
        self.time = time
        self.domain = domain
        self.is_1d = is_1d
