"""
This file is part of the AdaptiveMasterMSM package.

"""
import os
import subprocess
import shlex

# AdaptiveMasterMSM
from ..system import system_lib

def check_gmx():
    """
    Checks the available Gromacs version.
    
    """
    output = subprocess.check_output("which gmx", shell=True)
    print ("%s"%output)

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
        if self.gro is None: self.gro = 'conf.gro'
        gro = self.gro
        cmd = 'gmx editconf -f %s -o %s -c -d %f -bt cubic' %(gro, out, d)
        print (cmd)
        os.system(cmd)
        self.gro = out 

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
        self.gro = out

    def ionize(self, neutral=True, conc=None, inp='conf_solv.gro',\
                out='conf_solv_ions.gro', top='topol.top', tpr='ions.tpr'):
        """
        Adds ions for electroneutrality and / or target concentration.
        
        Parameters
        ----------
        neutral : bool
            Whether we neutralize the simulation box.
        conc : float
            The target salt concentration in M.

        """
        # format the water string
        if not conc:
            ion_str = '-neutral'
        else:
            ion_str = '-conc %g'%conc

        # generate the tpr file required for gmx genion
        self.mdpfile = system_lib.write_minimization_mdp(self.cons)
        cmd = "gmx grompp -f %s -c %s -p %s -o %s"%(self.mdpfile, inp, top, tpr)
        print(cmd)
        p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
        output, error = p.communicate()

        # runs genion 
        cmd = "gmx genion -s %s -o %s -p %s %s"%(tpr, out, top, ion_str)
        p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
        output, error = p.communicate()
        print(cmd)
        #return output, error
        self.gro = out

#    def minimize(self, inp='conf_solv_ions.gro', top='topol.top', out='minimization'):
#        """
#        Runs minimization invoking the Launcher class
#
#        """
#        self.mdpfile = system_lib.write_minimization_mdp(self.cons)
#        cmd = 'gmx grompp -f %s -c %s -p %s -o %s'% (self.mdpfile, inp, top, out)
#        print(cmd)
#        p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
#        output, error = p.communicate()
#
#        cmd = 'gmx mdrun -v -s %s.tpr -deffnm %s'%(out, out)
#        p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
#        output, error = p.communicate()
#        self.gro = out + ".gro"
#
#    def nvt(self, mdp='nvt.mdp', gro='minimization.gro',\
#            top='topol.top', tpr='nvt.tpr', out='nvt'):
#
#        cmd = 'gmx grompp -f %s -c %s -r %s -p %s -o %s; gmx mdrun -v -s %s -deffnm %s'\
#                % (mdp, gro, gro, top, tpr, tpr, out)
#        print(cmd)
#        p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
#        output, error = p.communicate()
#        return output, error
#
#    def npt(self, mdp='npt.mdp', gro='nvt.gro', chk='nvt.cpt',\
#            top='topol.top', tpr='npt.tpr', out='npt'):
#
#        cmd = 'gmx grompp -f %s -c %s -r %s -t %s -p %s -o %s; gmx mdrun -v -s %s -deffnm %s'\
#                % (mdp, gro, gro, chk, top, tpr, tpr, out)
#        print(cmd)
#        p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
#        output, error = p.communicate()
#        return output, error
