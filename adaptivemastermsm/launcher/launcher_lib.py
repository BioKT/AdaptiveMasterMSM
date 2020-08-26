"""
This file is part of the AdaptiveMasterMSM package.

"""
#!/usr/bin/env python

import sys, os

def checkfile(inp):

    if not os.path.exists(inp[:inp.rfind("/")]):
        os.makedirs(inp[:inp.rfind("/")])

def write_mdp(params):

    txt = """; GROMACS mdp options
    ; Run parameters
    integrator               = md
    dt                       = %f
    ; Output control
    nsteps                   = %d
    nstxout                  = 0
    nstvout                  = 0
    nstlog                   = %d
    nstenergy                = 0
    nstxout-compressed       = %d
    compressed-x-grps        = System
    ; Neighbors
    cutoff-scheme            = Verlet ; Buffered neighbor searching
    ns_type                  = grid   ; search neighboring grid cells
    rvdw                     = 0.1
    ; rvdw_switch              = 1.0
    pbc                      = xyz
    rlist                    = 0.03
    nstlist                  = 10
    ; Electrostatics
    coulombtype              = PME  ; Particle Mesh Ewald for long-range electrostatics
    rcoulomb                 = 0.1  ; short-range electrostatic cutoff (in nm)
    ; fourier_nx               = 0
    ; fourier_ny               = 0
    ; fourier_nz               = 0
    ; optimize_fft             = yes
    pme_order                = 4        ; cubic interpolation
    fourierspacing           = 0.08     ; grid spacing for FFT
    DispCorr                 = EnerPres ; account for cut-off vdW scheme
    ; T and P coupling
    tcoupl                   = %s
    tc_grps                  = System
    ; tc_grps                  = Protein Non-Protein ; two coupling groups-more accurate
    tau_t                    = 0.1
    ref_t                    = %f
    pcoupl                   = %s
    tau_p                    = 1
    ref_p                    = %f
    pcoupltype               = isotropic
    compressibility          = 0.000045
    ; Bond parameters
    gen_vel                  = no
    gen_temp                 = %f
    constraint_algorithm     = lincs    ; holonomic constraints
    constraints              = h-bonds  ; bonds involving H are constrained
    lincs_iter               = 1        ; accuracy of LINCS
    lincs_order              = 4        ; also related to accoutput
    continuation             = yes
    morse                    = no
    implicit_solvent         = no
    """ % ( params['timestep'],
    params['total_steps'],
    params['log_freq'],
    params['xtc_freq'],
    params['thermostat'],
    params['temperature'],
    params['barostat'],
    params['pressure'],
    params['temperature'] )

    return txt

def clean_working_directory():

    print ("Cleaning up...")
    os.system('rm \#*')
    os.system('mkdir logs')
    os.system('mv *.log logs')
    os.system('rm *.trr *.top *.itp *.edr *.cpt *.tpr out.gro conf.gro')
    os.system('mv equilibration.xtc *.mdp *.gro logs')
