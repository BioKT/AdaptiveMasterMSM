.. _documentation:

Modules
=============
Adaptivemastermsm is a Python package that is divided in three main subpackages. 
This way of structuring the code derives from the three main types of 
objects that are constructed. First, there is system, which 
result in objects of the ``System`` class; second, command execution
is managed by launcher, which come in the form of instances of the ``Launcher``
class; finally, MD data is analyzed by another objects, which are contained
in the ``Analyzer`` class.

SYSTEM module
-------------
This module contains everything necessary to get your time series data
into MasterMSM. The main class object within this module is the TimeSeries
object.

.. currentmodule:: adaptivemastermsm

.. autosummary::
    :toctree: 

    system


LAUNCHER module
---------------
.. currentmodule:: adaptivemastermsm

.. autosummary::
    :toctree:

    launcher


ANALYZER module
---------------
.. currentmodule:: adaptivemastermsm

.. autosummary::
    :toctree:

    analyzer
    
Examples
--------
We have put together a few simple Python notebooks to help you learn the basics
of the AdaptiveMasterMSM package. They are based on data derived from either model
systems or from molecular dynamics simulations of some simple (albeit realistic)
biomolecules. You can find the notebooks in the following 
`link <https://github.com/BioKT/AdaptiveMasterMSM/tree/develop/examples>`_.
