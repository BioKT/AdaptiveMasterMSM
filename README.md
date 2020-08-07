AdaptiveMasterMSM
=================
AdaptiveMasterMSM is a Python package for doing adaptive sampling by employing
Markov State Models (MSM) interacting with MasterMSM software and molecular
dynamics simulations through GROMACS. This package will allow you to:

*   Do adaptive sampling simulations

*   Employ different scoring functions: 'counts' and 'probabilities'

*   User MasterMSM for the MSM analysis, based on the master equation

You can read the documentation [here](https://adaptivemastermsm.readthedocs.io).

Contributors
------------
This code has been written by David De Sancho and Ion Mitxelena.

Installation
------------
    git clone http://github.com/BioKT/AdaptiveMasterMSM destination/AdaptiveMasterMSM
    cd destination/AdaptiveMasterMSM
    python setup.py install --user
    for developers use: python setup.py develop --user 

External libraries
------------------
    mdtraj : https://mdtraj.org
    MasterMSM: https://github.com/BioKT/MasterMSM
    GROMACS: http://www.gromacs.org
