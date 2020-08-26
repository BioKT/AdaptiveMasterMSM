"""
This file is part of the AdaptiveMasterMSM package.

"""
#!/usr/bin/env python

import sys, os
import subprocess
import shlex
import multiprocessing as mp

# AdaptiveMasterMSM
from adaptivemastermsm.system import system
from ..launcher import launcher_lib

class Launcher(object):
    """
    Call GROMACS and process output to pass it to class analyzer
    """

    def __init__(self, md_step, pdb, forcefield, water):
        """
        Read input from system, run MD, and process output
        for the class analyzer

        Parameters
        ----------
        md_step : str
            'Equilibration' or 'Production'
        pdb : str
            Path to the PDB or GRO file to start setup
        forcefield : str
            Name of force field in Gromacs format (e.g. amber03, charmm27...)
        water : str
            Water model for in Gromacs format (e.g. tip3p...)

        """

        self.pdb = pdb
        self.ff = forcefield
        self.wat = water

        if md_step == 'Equilibration':
            print ("md_step %s not valid" % md_step)
            raise Exception("Launcher only works with 'Production'")

    def runner(self, gmxinput):
        """
        Call MDRUN via multiprocessing

        """
        # set multiprocessing options
        n_threads = mp.cpu_count()
        pool = mp.Pool(processes=n_threads)
        # run simulations
        results = []
        for x in gmxinput:
            results.append(pool.apply_async(self.gromacs_worker, [x]))
        # close the pool and wait for each running task to complete
        pool.close()
        pool.join()
        for result in results:
            out, err = result.get()
            print("out: {} err: {}".format(out, err))

        launcher_lib.clean_working_directory()

    def inputs(self, gro_all, top_all=None, chk=None, mdp_params=None, mdp=None):
        """
        Prepare directories and files for all copies and launch all of them
        
        Parameters
        ----------
        gro_all : str
            List of paths to gro files corresponding to different
            resamplings separated by spaces only (omit commas)
        top_all : str
            List of paths to topology (top) files
        chk : str
            Path to checkpoint file
        mdp_params : dict
            Parameters for the MDP file defined on class Controller
        mdp_all : str
            Path to mdp file, only compatible with mdp_params=None

        """

        """
        !!! something here to filter bad or missing gro files,
        maybe this work if a single one is passed
        if len(gro_all.split()) == 1:
            gro = "data/gro/%s.gro" % i
            launcher_lib.checkfile(gro)
            cmd = "cp %s data/gro/%s.gro" % (self.pdb, i)
            os.system(cmd)
        """

        if top_all is None:
            top_all_list = []
            for i in range(len(gro_all.split())):
                top = "data/top/%s.top" % i
                launcher_lib.checkfile(top)
                top_all_list.append(top)
        else:
            top_all_list = top_all.split()

        # ndx ? use this to set positions of ligands in future work
        gmxinput = []
        i = 0
        for gro in gro_all.split():

            gro_out = "data/gro/processed/%s.gro" % i
            launcher_lib.checkfile(gro_out)
            top = top_all_list[i]
            top_out = "data/top/processed/%s.top" % i
            launcher_lib.checkfile(top_out)
            cmd = "cp %s %s" % (top, top_out)
            os.system(cmd)
            if top_all is None:
                cmd = 'gmx pdb2gmx -f %s -o %s -p %s -ff %s -water %s'%\
                      (gro, gro_out, top_out, self.ff, self.wat)
                print(cmd)
                os.system(cmd)
                itp = "data/itp/%s.itp" % i
                launcher_lib.checkfile(itp)
                cmd = "mv posre.itp %s" % itp
                os.system(cmd)
            else:
                cmd = "cp %s %s" % (gro, gro_out)
                os.system(cmd)
            if mdp is None: mdp = self.gen_mdp(i, mdp_params)
            out = "data/tpr/%s.tpr" % i
            self.gen_tpr(gro_out, top_out, mdp, out, i, chk=chk)
            gmxinput.append([out,gro_out])
            i += 1

        return gmxinput
    
    def gen_tpr(self, gro, top, mdp, tpr, i, chk=None):
        """
        Function for generating tpr files

        """
        launcher_lib.checkfile(tpr)
        if chk is not None:
            cmd = "gmx grompp -c %s -p %s -f %s -t %s -o %s -maxwarn 1"\
                %(gro, top, mdp, chk, tpr)
        else:
            cmd = "gmx grompp -c %s -p %s -f %s -o %s -maxwarn 1"\
                %(gro, top, mdp, tpr)
        print(cmd)
        os.system(cmd)
        cmd = "mv mdout.mdp data/mdp/%s_out.mdp" % i
        os.system(cmd)
    
    def gen_mdp(self, i, mdp_params):
        """ 
        Generate Gromacs parameter input file (mdp)
        
        """
        txt = launcher_lib.write_mdp(mdp_params)
        filemdpout = "data/mdp/%s.mdp" % i
        launcher_lib.checkfile(filemdpout)
        out = open(filemdpout, "w")
        out.write(txt)
        out.close()
        """ keep for a while to consider other options:
        rawmdp = open(filemdp).readlines()
        for i in rawmdp:
            out.write("%s"%i)
        """
        return filemdpout
    
    def gromacs_worker(self, x):
        """
        Worker function for running Gromacs jobs

        -->Shall we use 'multi' flag from Gromacs?

        """
        #print(x)
        tpr = x[0]
        out = x[1]

        # run simulation
        cmd = "gmx mdrun -v -s %s -deffnm %s"%(tpr, out)
        print(cmd)
        p = subprocess.Popen(shlex.split(cmd), stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        out, err = p.communicate()
        return (out, err)

    def run_md(self, n_threads, params):
        """
        worker provisional para 'Equilibration'
        """

        # note for the moment wet_xtc_file, dry_xtc_file and last_wet_snapshot
        # overwrite each other, until I see their need

        cmd = 'gmx grompp -f %s -c processed_%s -p topol.top -maxwarn 1'\
                % (params.filemdp, self.pdb) &\
                'gmx mdrun -nt %d -s topol.tpr -x %s -c processed_%s -g prod.log' % \
                (n_threads, self.pdb, self.wet_xtc_file)
        print(cmd)
        p = subprocess.Popen(shlex.split(cmd), stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        out, err = p.communicate()

        # if running in production mode, trjconv to create a dry XTC
        if params.md_step == 'Production':
            assert self.dry_xtc_file is not None and self.last_wet_snapshot is not None
            cmd = 'echo 0 | gmx trjconv -f end.gro -s topol.tpr -o %s -pbc whole' % self.last_wet_snapshot
            print(cmd)
            p = subprocess.Popen(shlex.split(cmd), stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            output, error = p.communicate()
            #print(output, error)
            cmd = 'echo PROTEIN | gmx trjconv -f %s -s topol.tpr -o %s -pbc whole' % (self.wet_xtc_file, self.dry_xtc_file)
            print(cmd)
            p = subprocess.Popen(shlex.split(cmd), stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            output, error = p.communicate()
            #print(output, error)

        return (out, err)
