.. hm documentation master file, created by
   sphinx-quickstart on Tue May 12 20:23:23 2020.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

hm: a Python framework for hydrological modelling
=================================================

**hm** is a free and open source project for building environmental models
in Python.

hm relies heavily on xarray_ for handling spatio-temporal datasets. The
underlying data model is provided by the netCDF_ file format. Methods
for timestepping are inspired by those included in PCRaster_.

.. _xarray: http://xarray.pydata.org
.. _netCDF: http://www.unidata.ucar.edu/software/netcdf
.. _PCRaster: http://pcraster.geo.uu.nl

Documentation
-------------

**Getting Started**

* :doc:`overview`
* :doc:`installing`

.. toctree::
   :maxdepth: 1
   :hidden:
   :caption: Getting Started

   overview
   examples
   installing

**User Guide**

* :doc:`data-structures`
* :doc:`configuration`
* :doc:`reporting`

.. toctree::
   :maxdepth: 1
   :hidden:
   :caption: User Guide

   data-structures
   configuration
   reporting

**Examples**

* :doc:`examples`
   
**Help & Reference**

.. toctree::
   :maxdepth: 2
   :caption: Contents:


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
