##################
AdaptiveMasterMSM
##################

.. image:: Logo_Web.pdf
   :width: 400

AdaptiveMasterMSM is an open-source implementation of adaptive sampling algorithms exploiting markov state models based on the rate matrix.

The AdaptiveMasterMSM python package is designed to carry out adaptive sampling
unbiased simulations in order to cover slow biological time scales. Our
implementation uses Gromacs and multiprocessing to simulate parallel short
trajectories. Markov state models are employed between consecutive trajectories
to determine where best to sample from. The latter is decided according to
(1) low-counts and (2) populations metrics, despite more options are being
investigated. In contrast to other similar computer programs, AdaptiveMasterMSM
employs the MasterMSM python package to construct the Markov state models,
which works with the rate matrix instead of the usually employed transition
matrix.

.. toctree::
   :maxdepth: 3
   :caption: Overview:
   
   Installation
   License
   Contact

.. toctree::
   :maxdepth: 3
   :caption: User Documentation:
   
   Getting-Started

.. Indices and tables
.. ==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

