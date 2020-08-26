"""
This file is part of the AdaptiveMasterMSM package.

"""
#!/usr/bin/env python

import os, sys

def write_minimization_mdp(constraints='None'):

    txt = """; GROMACS mdp minimization options
    integrator               = steep
    emstep                   = 0.001
    emtol                    = 10.0
    nsteps                   = 50000
    nstxout                  = 0
    nstvout                  = 0
    nstlog                   = 10
    nstenergy                = 0
    nstxtcout                = 0
    xtc_grps                 = System
    energygrps               = 
    nstlist                  = 1
    ns_type                  = grid
    pbc                      = xyz
    periodic_molecules       = no
    rlist                    = 1.0
    rcoulomb                 = 1.0
    rvdw                     = 1.0
    tcoupl                   = no
    Pcoupl                   = no
    gen_vel                  = yes
    constraints              = %s
    continuation             = no
    morse                    = no
    implicit_solvent         = no
    """ % constraints
        
    filemdp = "data/minimization.mdp"
    try:
        f = open(filemdp, 'w')
        f.write(txt)
    except IOError:
        os.makedirs(filemdp[:filemdp.rfind("/")])
        f = open(filemdp, "w")
        f.write(txt)
    f.close()

    return filemdp
