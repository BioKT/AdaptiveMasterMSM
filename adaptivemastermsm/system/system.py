"""
This file is part of the AdaptiveMasterMSM package.

"""
#!/usr/bin/env python

import os
import subprocess
import shlex

# AdaptiveMasterMSM
from ..system import system_lib

class System(object):
    """
    Create files to be able to run GROMACS
    """

    def __init__(self, pdb, forcefield=None, water=None,\
                 topology=None, constraints=None):
        """
        Parameters
        ----------
        pdb : str
            Path to the PDB or GRO file to start setup
        forcefield : str
            Name of force field in Gromacs format (e.g. amber96)
        water : str
            Water model for GROMACS (e.g. tip3p)
        topology : str
            Path to topology TOP file
        constraints : str
            Parameter for the minimization MDP file (e.g. h-bonds)

        """
        # Read user-defined parameters
        self.pdb = pdb
        if forcefield is not None:
            self.ff = forcefield
            self.wat = water
        elif topology is not None:
            self.top = topology
        self.cons = constraints
    
    def gen_top(self, gro='conf.gro', top='topol.top'):
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

    def solvate(self, inp='conf_edit.gro', out='conf_solv.gro', top='topol.top'):
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

        cmd = 'gmx solvate -cp %s -cs %s -p %s -o %s' %\
                (inp, water_topol, top, out)
        print (cmd)
        os.system(cmd)

    def ionize(self, anion=0, cation=0, inp='conf_solv.gro',\
                out='conf_solv_ions.gro', top='topol.top', tpr='ions.tpr'):
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
        self.mdpfile = system_lib.write_minimization_mdp(self.cons)
        cmd = "gmx grompp -f %s -c %s -p %s -o %s; gmx genion -s %s -o %s -p %s %s"\
                % (self.mdpfile, inp, top, tpr, tpr, out, top, ion_str)
        print(cmd)
        p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
        output, error = p.communicate()
        return output, error
 
    def minimize(self, inp='conf_solv_ions.gro', top='topol.top'):

        self.mdpfile = system_lib.write_minimization_mdp(self.cons)
        cmd = 'gmx grompp -f %s -c %s -p %s -o minimization; gmx mdrun -v -s minimization -x minimization.xtc -deffnm minimization'\
                % (self.mdpfile, inp, top)
        print(cmd)
        p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
        output, error = p.communicate()
        return output, error

    def nvt(self, mdp='nvt.mdp', gro='minimization.gro',\
            top='topol.top', tpr='nvt.tpr', out='nvt'):

        cmd = 'gmx grompp -f %s -c %s -r %s -p %s -o %s; gmx mdrun -v -s %s -deffnm %s'\
                % (mdp, gro, gro, top, tpr, tpr, out)
        print(cmd)
        p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
        output, error = p.communicate()
        return output, error

    def npt(self, mdp='npt.mdp', gro='nvt.gro', chk='nvt.cpt',\
            top='topol.top', tpr='npt.tpr', out='npt'):

        cmd = 'gmx grompp -f %s -c %s -r %s -t %s -p %s -o %s; gmx mdrun -v -s %s -deffnm %s'\
                % (mdp, gro, gro, chk, top, tpr, tpr, out)
        print(cmd)
        p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
        output, error = p.communicate()
        return output, error
