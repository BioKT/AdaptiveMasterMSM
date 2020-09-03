"""
This file is part of the AdaptiveMasterMSM package.

"""
#!/usr/bin/env python

import sys, os

def checkfile(inp):

    if not os.path.exists(inp[:inp.rfind("/")]):
        os.makedirs(inp[:inp.rfind("/")])

def write_mdp_min():

    txt = """
    ; GROMACS mdp minimization options
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
    constraints              = none
    continuation             = no
    morse                    = no
    implicit_solvent         = no
    """
    return txt

def write_mdp_nvt():

    txt = """
    define                  = -DPOSRES
    integrator              = md
    dt                      = 0.002 ; ps
    nsteps                  = 50000 ; 100 ps
    nstlog                  = 1000
    nstenergy               = 1000
    nstxout-compressed      = 1000
    constraint_algorithm    = lincs
    constraints             = h-bonds
    lincs_iter              = 1
    lincs_order             = 4
    nstlist                 = 10
    ns_type                 = grid
    cutoff-scheme           = verlet
    coulombtype             = pme
    rvdw                    = 1.
    rcoulomb                = 1.
    ;# with PME we should not have to use coupling groups
    tcoupl                  = v-rescale
    tc_grps                 = protein non-protein
    tau_t                   = 0.1   0.1
    ref_t                   = 300   300
    Pcoupl                  = no
    pbc                     = xyz
    gen_vel                 = yes
    gen_temp                = 10.0
    gen_seed                = -1
    """
    return txt

def write_mdp_npt():

    txt = """
    refcoord_scaling        = com
    integrator              = md
    dt                      = 0.002
    nsteps                  = 50000
    nstlog                  = 100
    nstenergy               = 100
    nstxout-compressed      = 100
    compressed-x-grps       = non-Water
    continuation            = yes
    constraint_algorithm    = lincs
    constraints             = h-bonds
    lincs_iter              = 1
    lincs_order             = 4
    cutoff-scheme           = verlet
    nstlist                 = 10
    ns_type                 = grid
    coulombtype             = pme
    pme_order               = 4
    fourier_spacing         = 0.16
    rvdw                    = 1.
    rcoulomb                = 1.
    tcoupl                  = v-rescale
    tc_grps                 = protein non-protein
    tau_t                   = 0.1   0.1
    ref_t                   = 300   300
    ; Isotropic pressure coupling is now on
    Pcoupl                  = Berendsen
    Pcoupltype              = isotropic
    tau_p                   = 2.0
    ref_p                   = 1.0
    compressibility         = 4.5e-5
    pbc                     = xyz
    gen_vel                 = no
    """
    return txt

def write_mdp_user(params):

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
