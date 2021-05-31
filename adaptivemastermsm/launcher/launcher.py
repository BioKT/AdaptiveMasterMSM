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

    def gen_tpr(self, mdp=None, tpr='data/out.tpr', nsteps=250000):#, sim_time=None,chk=None):
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
        emstep = 0.002
        if mdp in ["min","nvt","npt","prod"]:
            if mdp is "min": txt = launcher_lib.write_mdp_min()
            if mdp is "nvt": txt = launcher_lib.write_mdp_nvt()
            if mdp is "npt": txt = launcher_lib.write_mdp_npt()
            if mdp is "prod": txt = launcher_lib.write_mdp_prod(emstep=emstep, nsteps=nsteps)
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
        #sim_time += emstep * float(nsteps)

    def gmx_run(self, out='data/out'):
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

    def mp_run(self, tprs, outs):
        """ Driver for gmx_worker

        """
        # define multiprocessing options
        nproc = mp.cpu_count()
        if len(tprs) > nproc:
            sys.exit(" Please reduce number of parallel runs")
        #self.nt = nproc/len(tprs)

        pool = mp.Pool(processes=nproc)
        # generate multiprocessing input
        mpinput = []
        for i, tpr in enumerate(tprs):
            mpinput.append([tpr, outs[i]])
        #mpinput = [[tpr, out] for tpr in tprs for out in outs]
        # run counting using multiprocessing
        pool.map(self.gmx_worker, mpinput)
        pool.close()
        pool.join()

    def gmx_worker(self, x):
        """ mp worker for running Gromacs jobs

        -->Shall we use 'multi' flag from Gromacs?

        Parameters
        ----------
        x : list
            List containing input for each mp worker. Includes:
            path to tpr file
            path to inp/out file

        """
        #print(x)
        tpr = x[0]
        out = x[1]

        # run simulation
        #if self.nt > 1:
        #cmd = "gmx mdrun -nt %g -v -s %s -deffnm %s"%(self.nt, tpr, out)
        #else:
        cmd = "gmx mdrun -nt 1 -v -s %s -deffnm %s"%(tpr, out)
        print(cmd)
        p = subprocess.Popen(shlex.split(cmd), stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        out, err = p.communicate()
        #return (out, err)