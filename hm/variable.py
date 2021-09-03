#!/usr/bin/env python
# -*- coding: utf-8 -*-

from typing import List
# See the following link for background on typing.List vs built-in list
# https://stackoverflow.com/a/39458225
# (tl;dr - typing.List allows the user to specify types of contents)

class Variable:
    """Represents a model input/output variable and associated metadata."""

    def __init__(
            self,
            standard_name: str,
            variable: str,
            variable_type: str,
            variable_item_size: int,
            variable_class: str,
            units: str,
            grid: int
    ):
        self.standard_name: str = standard_name
        self.variable: str = variable
        self.variable_type: str = variable_type
        self.variable_item_size: int = variable_item_size
        self.variable_class: str = variable_class
        self.units: str = units
        self.grid: int = grid

def get_variable_list(variable_dict):
    variable_list: List[Variable] = []
    for var in variable_dict:
        variable_list.append(
            Variable(
                standard_name = var['standard_name'],
                variable = var['variable'],
                variable_type = var['variable_type'],
                variable_item = var['variable_item_size'],
                variable_class = var['variable_class'],
                units = var['units'],
                grid = var['grid']
            )
        )
            
    return variable_list
        
class IO(object):
    def __init__(
            self,
            input_variable_dict,
            output_variable_dict# ,
            # state_variable_dict
    ):
        self.inputs = get_variable_list(input_variable_dict)
        self.output = get_variable_list(output_variable_dict)
        # self.states = get_variable_list(state_variable_dict)

    @property
    def inputs(self):
        return self.inputs

    @property
    def outputs(self):
        return self.outputs

    @property
    def inputs_outputs(self):
        return self.inputs + self.outputs
        
