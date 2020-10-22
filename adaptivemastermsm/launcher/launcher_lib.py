"""
This file is part of the AdaptiveMasterMSM package.

"""
#!/usr/bin/env python

import sys, os

def checkfile(inp):

    if not os.path.exists(inp[:inp.rfind("/")]):
        os.makedirs(inp[:inp.rfind("/")])

#def clean_working_directory():
#
#    print ("Cleaning up...")
#    os.system('rm \#*')
#    os.system('mkdir logs')
#    os.system('mv *.log logs')
#    os.system('rm *.trr *.top *.itp *.edr *.cpt *.tpr out.gro conf.gro')
#    os.system('mv equilibration.xtc *.mdp *.gro logs')

def write_mdp_min(emstep=0.001, nsteps=50000, constraints="none"):

    txt = """
    ; GROMACS mdp minimization options
    integrator               = steep
    emstep                   = %f
    emtol                    = 10.0
    nsteps                   = %d
    nstxout                  = 0
    nstvout                  = 0
    nstlog                   = 1000
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
    """ % (emstep, nsteps, constraints)
    return txt

def write_mdp_nvt(dt=0.002, nsteps=50000, constraints="h-bonds", ref_t=300):

    txt = """
    define                  = -DPOSRES
    integrator              = md
    dt                      = %f ; ps
    nsteps                  = %d ; 100 ps
    nstlog                  = 1000
    nstenergy               = 1000
    nstxout-compressed      = 1000
    constraint_algorithm    = lincs
    constraints             = %s
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
    ref_t                   = %d   %d
    Pcoupl                  = no
    pbc                     = xyz
    gen_vel                 = yes
    gen_temp                = 10.0
    gen_seed                = -1
    """ % (dt, nsteps, constraints, ref_t, ref_t)
    return txt

def write_mdp_npt(dt=0.002, nsteps=50000, constraints="h-bonds",\
                    ref_t=300, ref_p=1.0):

    txt = """
    refcoord_scaling        = com
    integrator              = md
    dt                      = %f
    nsteps                  = %d
    nstlog                  = 1000
    nstenergy               = 1000
    nstxout-compressed      = 1000
    compressed-x-grps       = non-Water
    continuation            = yes
    constraint_algorithm    = lincs
    constraints             = %s
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
    ref_t                   = %d   %d
    ; Isotropic pressure coupling is now on
    Pcoupl                  = Berendsen
    Pcoupltype              = isotropic
    tau_p                   = 2.0
    ref_p                   = %f
    compressibility         = 4.5e-5
    pbc                     = xyz
    gen_vel                 = no
    """ % (dt, nsteps, constraints, ref_t, ref_t, ref_p)
    return txt

def write_mdp_prod(emstep=0.002, nsteps=50000, constraints="h-bonds"):

    txt = """
    ; GROMACS mdp production options
    integrator   		    = sd
    dt                      = %f ; in ps
    nsteps  			    = %d
    nstlog			        = 100
    nstenergy       		= 100
    nstxout-compressed	    = 100
    ; nstxtcout               = 100
    ; xtc-grps                = protein
    compressed-x-grps	    = non-Water
    continuation		    = yes
    constraint_algorithm    = lincs  
    constraints             = %s
    lincs_iter              = 1      
    lincs_order             = 4      
    nstlist			        = 10
    ns_type			        = grid
    cutoff-scheme		    = verlet 
    coulombtype 		    = pme
    rvdw                    = 1.
    rcoulomb               	= 1.
    tc_grps                 = protein non-protein 
    tau_t                   = 1	1
    ref_t                   = 300.0 300.0
    Pcoupl                  = no
    gen_vel                 = no
    """ % (emstep, nsteps, constraints)
    return txt

#integrator              = md
#dt                      = 0.004
#nsteps                  = 150000000 ; 600ns
#nstlog                  = 2500
#nstxout                 = 2500
#nstvout                 = 2500
#nstfout                 = 2500
#nstcalcenergy           = 2500
#nstenergy               = 2500
#;energygrps              = Protein Non-Protein ; Not available for GPU-based runs
#;
#cutoff-scheme           = Verlet
#nstlist                 = 20
#coulombtype             = pme
#rcoulomb                = 1.2
#vdwtype                 = Cut-off
#vdw-modifier            = Force-switch
#rvdw_switch             = 1.0
#rvdw                    = 1.2
#;
#tcoupl                  = V-rescale
#tc_grps                 = Protein Non-Protein
#tau_t                   = 0.1  0.1
#ref_t                   = 310   310
#nsttcouple              = 10
#;
#pcoupl                  = Parrinello-Rahman
#pcoupltype              = isotropic
#tau_p                   = 2.0
#compressibility         = 4.5e-5
#ref_p                   = 1.0
#;
#constraints             = all-bonds
#constraint_algorithm    = LINCS
#continuation            = yes
#;
#refcoord_scaling        = com
#DispCorr                = EnerPres
#gen-vel                 = no     ; Assign velocities to particles by taking them randomly from a Maxwell distributioni
#gen-temp                = 310

#def write_mdp_user(params):
#
#    txt = """; GROMACS mdp options
#    ; Run parameters
#    integrator               = md
#    dt                       = %f
#    ; Output control
#    nsteps                   = %d
#    nstxout                  = 0
#    nstvout                  = 0
#    nstlog                   = %d
#    nstenergy                = 0
#    nstxout-compressed       = %d
#    compressed-x-grps        = System
#    ; Neighbors
#    cutoff-scheme            = Verlet ; Buffered neighbor searching
#    ns_type                  = grid   ; search neighboring grid cells
#    rvdw                     = 0.1
#    ; rvdw_switch              = 1.0
#    pbc                      = xyz
#    rlist                    = 0.03
#    nstlist                  = 10
#    ; Electrostatics
#    coulombtype              = PME  ; Particle Mesh Ewald for long-range electrostatics
#    rcoulomb                 = 0.1  ; short-range electrostatic cutoff (in nm)
#    ; fourier_nx               = 0
#    ; fourier_ny               = 0
#    ; fourier_nz               = 0
#    ; optimize_fft             = yes
#    pme_order                = 4        ; cubic interpolation
#    fourierspacing           = 0.08     ; grid spacing for FFT
#    DispCorr                 = EnerPres ; account for cut-off vdW scheme
#    ; T and P coupling
#    tcoupl                   = %s
#    tc_grps                  = System
#    ; tc_grps                  = Protein Non-Protein ; two coupling groups-more accurate
#    tau_t                    = 0.1
#    ref_t                    = %f
#    pcoupl                   = %s
#    tau_p                    = 1
#    ref_p                    = %f
#    pcoupltype               = isotropic
#    compressibility          = 0.000045
#    ; Bond parameters
#    gen_vel                  = no
#    gen_temp                 = %f
#    constraint_algorithm     = lincs    ; holonomic constraints
#    constraints              = h-bonds  ; bonds involving H are constrained
#    lincs_iter               = 1        ; accuracy of LINCS
#    lincs_order              = 4        ; also related to accoutput
#    continuation             = yes
#    morse                    = no
#    implicit_solvent         = no
#    """ % ( params['timestep'],
#    params['total_steps'],
#    params['log_freq'],
#    params['xtc_freq'],
#    params['thermostat'],
#    params['temperature'],
#    params['barostat'],
#    params['pressure'],
#    params['temperature'] )
#
#    return txt
#

