"""
This file is part of the AdaptiveMasterMSM package.

"""
#!/usr/bin/env python

import sys, os
import subprocess
import shlex
import multiprocessing as mp

def checkfile(inp):

    if not os.path.exists(inp[:inp.rfind("/")]):
        os.makedirs(inp[:inp.rfind("/")])

    return

def write_mdp():

    txt = """; GROMACS mdp options
    integrator               = md
    dt                       = %f
    nsteps                   = %d
    nstxout                  = 0
    nstvout                  = 0
    nstlog                   = %d
    nstenergy                = 0
    nstxtcout                = %d
    xtc_grps                 = System
    nstlist                  = 1
    ns_type                  = grid
    pbc                      = xyz
    periodic_molecules       = no
    rlist                    = 0.03
    coulombtype              = PME
    fourier_nx               = 0
    fourier_ny               = 0
    fourier_nz               = 0
    optimize_fft             = yes
    pme_order                = 4
    fourierspacing           = 0.08
    rcoulomb                 = 0.1
    vdwtype                  = shift
    rvdw                     = 0.1
    ; rvdw_switch              = 1.0
    tcoupl                   = %s
    tc_grps                  = System
    tau_t                    = 1
    ref_t                    = %f
    Pcoupl                   = %s
    tau_p                    = 1
    pcoupltype               = isotropic
    compressibility          = 0.000045
    ref_p                    = %f
    gen_vel                  = no
    gen_temp                 = %f
    constraints              = hbonds
    continuation             = yes
    ; continuation             = no
    morse                    = no
    implicit_solvent         = no
    """ % ( self.run['timestep'],
    self.run['total_steps'],
    self.run['log_freq'],
    self.run['xtc_freq'],
    self.run['thermostat'],
    self.run['temperature'],
    self.run['barostat'],
    self.run['pressure'],
    self.run['temperature'] )

    return txt
