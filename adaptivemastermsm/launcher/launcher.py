"""
This file is part of the AdaptiveMasterMSM package.

"""
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

    def __init__(self, simsys):
        """
        Read input from system, run MD, and process output
        for the class analyzer

        Parameters
        ----------
        simsys : object
            An instance of the System class.

        """
        Launcher.system = simsys

    def gen_tpr(self, mdp=None, tpr='out.tpr'): #, i, chk=None):
        """
        Function for generating tpr files

        Parameters
        ----------
        mdp : str
            The parameter input file
        tpr : str
            The Gromacs run input file

        """
        #launcher_lib.checkfile(tpr)
        #if chk is not None:
        #    cmd = "gmx grompp -c %s -p %s -f %s -t %s -o %s -maxwarn 1"\
        #        %(gro, top, mdp, chk, tpr)
        #else:
        filemdp = mdp
        if mdp in ["min","nvt","npt"]:
            if mdp is "min": txt = launcher_lib.write_mdp_min()
            if mdp is "nvt": txt = launcher_lib.write_mdp_nvt()
            if mdp is "npt": txt = launcher_lib.write_mdp_npt()
            filemdp = "data/mdp/%s.mdp" % mdp
            launcher_lib.checkfile(filemdp)
            inp = open(filemdp, "w")
            inp.write(txt)
            inp.close()
        cmd = "gmx grompp -c %s -p %s -f %s -o %s -r %s"\
                %(self.system.gro, self.system.top, filemdp, tpr, self.system.gro)
        print(cmd)
        os.system(cmd)
        self.tpr = tpr

    def gmx_run(self, out='out'):
        """
        Runs MD simulations using gmx mdrun

        Parameters
        ----------
        tpr : str
            Gromacs run input file
        out : str
            Root filename for output

        """
        cmd = "gmx mdrun -v -s %s -deffnm %s"%(self.tpr, out)
        print(cmd)
        p = subprocess.Popen(shlex.split(cmd), \
                stdout=subprocess.PIPE, \
                stderr=subprocess.PIPE)
        output, error = p.communicate()

        self.output = output
        self.error = error
        self.out = out
        # after the run we replace the gro in the system for the output gro file
        self.system.gro =  "%s.gro"%out

#    def runner(self, tpr, out):
#        """
#        Call MDRUN via multiprocessing
#
#        """
#        # set multiprocessing options
#        n_threads = mp.cpu_count()
#        pool = mp.Pool(processes=n_threads)
#        # run simulations
#        results = []
#        for x in gmxinput:
#            results.append(pool.apply_async(self.gromacs_worker, [x]))
#        # close the pool and wait for each running task to complete
#        pool.close()
#        pool.join()
#        for result in results:
#            out, err = result.get()
#            print("out: {} err: {}".format(out, err))
#
#        launcher_lib.clean_working_directory()

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
    
    def gen_mdp(self, i, mdp_params):
        """ 
        Generate Gromacs parameter input file (mdp)
        
        """
        txt = launcher_lib.write_mdp_user(mdp_params)
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