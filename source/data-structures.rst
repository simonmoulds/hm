Data Structures
===============

`hm` includes three main data structures for representing data
in environmental models. 

HmDataArray
-----------

HmSpaceDataArray
^^^^^^^^^^^^^^^^

:py:class:`hm.HmSpaceDataArray` is a wrapper for `xarray.DataArray`
for spatial datasets which have no time dimension. Examples from
environmental modelling applications include maps of topography,
land use/land cover, and soil properties. 

HmSpaceTimeDataArray
^^^^^^^^^^^^^^^^^^^^

:py:class:`hm.SpaceTimeDataArray` is a wrapper for `xarray.DataArray`
for spatiotemporal datasets. This class is used to represent data
such as meteorological input data.

Creating a HmDataArray
----------------------

The constructor function `hm.open_hmdataarray` takes the following
arguments:

- ``filename_or_obj``
- ``variable``
- ``domain``
- ``is_1d``
- ``xy_dimname``
- ``use_xarray``
- ``xarray_kwargs``
- ``**kwargs``

HmModelTime
-----------

:py:class:`hm.ModelTime` represents the temporal grid of a
hydrological model.

Creating a HmModelTime

The constructor function `hm.open_modeltime` takes the following
arguments:

- ``starttime``
- ``endtime``
- ``timedelta``

HmDomain
--------

:py:class:`hm.Domain` is a class to represent the spatiotemporal
model domain.

Creating a HmModelDomain
------------------------

The constructor function `hm.set_domain` takes the following
arguments:

- ``filename_or_obj``
- ``modeltime``
- ``mask_varname``
- ``area_varname``
- ``is_1d``
- ``xy_dimname``
- ``pseudo_coords``
- ``**kwargs``
