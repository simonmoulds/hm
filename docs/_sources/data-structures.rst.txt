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

.. autofunction:: hm.api.open_hmdataarray

HmModelTime
-----------

:py:class:`hm.ModelTime` represents the temporal grid of a
hydrological model.

Creating a HmModelTime
----------------------

.. autofunction:: hm.api.open_modeltime

HmDomain
--------

.. autoclass:: hm.dataarray.HmSpaceDataArray

Creating a HmModelDomain
------------------------

.. autofunction:: hm.api.set_domain
		  
