"""
This file is part of the AdaptiveMasterMSM package.

"""
#!/usr/bin/env python

import os
import subprocess
import shlex

class System(object):
    """
    Create files to be able to run GROMACS
    """

    def __init__(self, pdb, forcefield=None, water=None, topology=None):
        """
        Parameters
        ----------
        pdb : str
            Path to the PDB or GRO file to start setup
        forcefield : str
            Name of force field in Gromacs format (e.g. amber03, charmm27...)
        water : str
            water model for GROMACS (tip3p, spc, ...)
        top : str
            Path to topology TOP file

        """
        # Read user-defined parameters
        self.pdb = pdb
        if forcefield is not None:
            self.ff = forcefield
            self.wat = water
        elif topology is not None:
            self.top = topology
    
    def gen_top(self, gro='conf', top='topol'):
        """
        Generates topology from structure file.

        """
        pdb = self.pdb
        self.gro = gro
        ff = self.ff
        wat = self.wat
        self.top = top
        cmd = "gmx pdb2gmx -f %s -o %s -p %s -ff %s -water %s"%\
                (pdb, gro, top, ff, wat)
        print (cmd)
        os.system(cmd)

    def build_box(self, out='conf_edit', d=1.0):
        """
        Generates the box for the system

        Parameters
        ----------
        d : float
            Distance between solute and box sides        
        
        """
        gro = self.gro
        cmd = 'gmx editconf -f %s -o %s -bt cubic -d %f -center' %(gro, out, d)
        print (cmd)
        os.system(cmd)

    def solvate(self):
        """
        Solvates the system using gmx solvate

        """
        # choose a water topology file
        water_dict  = {'tip3p' : 'spc216.gro'}
        if self.water in water_dict.keys():
            water_topol = water_dict[self.water]
        else:
            print ("Could not find a topology in 'w_top_dict' for water: %s" % self.water)
            raise Exception("Invalid water model for GROMACS.")

        'gmx genbox -cp out.gro -cs %s -p topol.top' % water_topol

#    def ionize(self):
#        """
#        Adds ions for electroneutrality
#
#        """
#        # format the water string
#        ion_str = ''
#        if self.run['Cl'] > 0:
#            ion_str += '-nn %d ' % self.run['Cl']
#        if self.run['Na'] > 0:
#            ion_str += '-np %d ' % self.run['Na']
 
#    def minimize(self)
#
#        mdpfile = self.write_minimization_mdp(i)
#
#        'gmx grompp -f %s -c out.gro -p topol.top' % mdpfile; \
#        'echo SOL | gmx genion -s topol.tpr -o out.gro -p topol.top %s' % ion_str; \
#        'gmx grompp -f %s -c out.gro -p topol.top' % mdpfile; \
#        'gmx mdrun -v -s topol.tpr -x minimization.xtc -c processed_%s -g EM.log' % filepdb
#        
#        print(" running: ",cmd)
#        p = subprocess.Popen(shlex.split(cmd), stdout=subprocess.PIPE, stderr=subprocess.PIPE)
#        out, err = p.communicate()
#        return out, err


    def runner():
        self.md_step = md_step

        # Take production or equilibration paths
        if self.md_step == 'Equilibration':
            self.driver_equilibration()
        elif self.md_step == 'Production':
            self.driver_production()
        else:
            print ("md_step %s not valid" % self.md_step)
            raise Exception("It must be 'Production' or 'Equilibration'")

        self.filemdp = self.write_mdp(i)

        return
        
    def driver_equilibration(self):
    
        self.equil = self.set_equilibration()
        self.run = {}
        self.equil['total_steps'] = \
            int(self.equil['total_time'] / self.equil['timestep'])
        
        for p in self.equil.keys():
            self.run[p] = self.equil[p]

        return           

    def driver_production(self):

        self.prod = self.set_production()
        self.run = {}
        self.prod['total_steps'] = \
            int(self.prod['total_time'] / self.prod['timestep'])

        for p in self.prod.keys():
            self.run[p] = self.prod[p]

        return

    def set_production(self):
        """
        Parameters for the MD production run
        """
        production = {
            'timestep'       : 0.002,         # ps
            'total_time'     : 1,             # ps / 10ns
            'log_freq'       : 10000,         # timesteps / 20ps
            'xtc_freq'       : 100,         # timesteps
            'temperature'    : 300,           # Kelvin
            'thermostat'     : 'Nose-Hoover', # type
            'box_size'       : 3.5,           # Angstroms
            'barostat'       : 'No',          # No barostat
            'pressure'       : 1,             # pressure
            'Cl'             : 0,             # number Cl's
            'Na'             : 0              # number Na's
        }

        return production

    def set_equilibration(self):
        """
        Parameters for the MD equilibration run
        """
        equilibration = {
            'timestep'       : 0.002,         # ps
            'total_time'     : 1,          # ps
            'barostat'       : 'berendsen',   # type
            'pressure'       : 1,              # bar
            'log_freq'       : 10000,         # timesteps / 20ps
            'xtc_freq'       : 100,         # timesteps
            'temperature'    : 300,           # Kelvin
            'thermostat'     : 'Nose-Hoover', # type
            'box_size'       : 3.5,           # Angstroms
            'Cl'             : 0,             # number Cl's
            'Na'             : 0              # number Na's
        }

        return equilibration

   
    def write_mdp(self, i):

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

        # ion: poner flag para check si existe ya el file
        filemdp = "data/as_%g.mdp" % i
        try:
            f = open(filemdp, 'w')
            f.write(txt)
        except IOError:
            os.makedirs(filemdp[:filemdp.rfind("/")])
            f = open(filemdp, "w")
            f.write(txt)
        f.close()

        return filemdp

    def write_minimization_mdp(self, i):

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
    constraints              = none
    continuation             = no
    morse                    = no
    implicit_solvent         = no
    """
        
        # ion: poner flag para check si existe ya el file
        filemdp = "data/as_minimization_%g.mdp" % i
        try:
            f = open(filemdp, 'w')
            f.write(txt)
        except IOError:
            os.makedirs(filemdp[:filemdp.rfind("/")])
            f = open(filemdp, "w")
            f.write(txt)
        f.close()

        return filemdp

