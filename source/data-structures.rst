Data Structures
===============

`hm` includes three main data structures for representing data
in environmental models. 

HmDataArray
-----------

HmSpaceDataArray
^^^^^^^^^^^^^^^^

.. autoclass:: hm.dataarray.HmSpaceDataArray

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

.. autoclass:: hm.dataarray.HmSpaceDataArray

Creating a HmModelDomain
------------------------

.. autofunction:: hm.api.set_domain
		  
