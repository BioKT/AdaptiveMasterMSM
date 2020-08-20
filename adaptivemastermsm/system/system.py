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
            Water model for GROMACS (tip3p, spc, ...)
        top : str
            Path to topology TOP file
        copy: int
            Number to identify which copy of the system we are dealing with

        """
        # Read user-defined parameters
        self.pdb = pdb
        if forcefield is not None:
            self.ff = forcefield
            self.wat = water
        elif topology is not None:
            self.top = topology
        #self.copy = copy
    
    def gen_top(self, gro='conf.gro', top='topol'):
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

    def build_box(self, out='conf_edit.gro', d=1.0):
        """
        Generates the box for the system

        Parameters
        ----------
        d : float
            Distance between solute and box sides        
        
        """
        gro = self.gro
        cmd = 'gmx editconf -f %s -o %s -c -d %f -bt cubic' %(gro, out, d)
        print (cmd)
        os.system(cmd)

    def solvate(self, inp='conf_edit.gro', out='conf_solv.gro'):
        """
        Solvates the system using gmx solvate

        """
        # choose a water topology file
        water_dict  = {'tip3p' : 'spc216.gro'}
        if self.wat in water_dict.keys():
            water_topol = water_dict[self.wat]
        else:
            print ("Could not find a topology in 'w_top_dict' for water: %s" % self.wat)
            raise Exception("Invalid water model for GROMACS.")

        cmd = 'gmx solvate -cp %s -cs %s -p topol.top -o %s' % (inp, water_topol, out)
        print (cmd)
        os.system(cmd)

    def ionize(self, anion=0, cation=0, inp='conf_solv.gro', out='conf_solv_ions.gro'):
        """
        Adds ions for electroneutrality
        
        Parameters
        ----------
        anion: int
            self.run['Cl']
        cation: int
            self.run['Na']

        """
        # format the water string
        ion_str = ''
        if anion > 0:
            ion_str += '-nn %d ' % anion
        if cation > 0:
            ion_str += '-np %d ' % cation
        mdpfile = self.write_minimization_mdp()
        cmd = 'gmx grompp -f %s -c %s -p topol.top -o ions.tpr' % (mdpfile, inp) &\
        'gmx genion -s ions.tpr -o %s -p topol.top %s' % (out, ion_str)
        print(" running: ",cmd)
        p = subprocess.Popen(shlex.split(cmd), stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        output, error = p.communicate()
        return output, error
 
    def minimize(self, inp='conf_solv_ions.gro'):

        mdpfile = self.write_minimization_mdp()
        cmd = 'gmx grompp -f %s -c %s -p topol.top -o minimization & gmx mdrun -v -s minimization -x minimization.xtc -deffnm minimization'\
                % (mdpfile, inp)
        print(" running: ",cmd)
        p = subprocess.Popen(shlex.split(cmd), stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        output, error = p.communicate()
        return output, error

    def write_minimization_mdp(self):

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
