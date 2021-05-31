"""
This file is part of the AdaptiveMasterMSM package.

"""
# Tools
import copy
import h5py
# Maths
import numpy as np
# Plotting
import matplotlib.pyplot as plt
# Multiprocessing
import multiprocessing as mp
# AdaptiveMasterMSM
from adaptivemastermsm import controller_lib
from adaptivemastermsm.system import system
from adaptivemastermsm.launcher import launcher
from adaptivemastermsm.launcher import launcher_lib
from adaptivemastermsm.analyzer import analyzer

class Controller(object):
    """
    Driver of the adaptive sampling algorithm

    """
    def __init__(self, pdb=None, forcefield=None, water=None,\
                topology=None, constraints=None, trajfiles=None, tprs=None):
        """
        Iteratively, create a GROMACS file, launch, and analyze.

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
        trajfiles : list
            Path(s) to trajectory file(s)
        tprs : list
            Path(s) to tpr file(s) corresponding to trajfiles

        """

        # Define the system
        self.gro_initial = pdb
        Controller.system = system.System(pdb, forcefield, water,\
                            topology=topology, constraints=constraints)

        self.top = topology
        if trajfiles is None:
            self.sys_prepare(pdb)
            self.all_sys_equilibration()
            self.trajfiles = ['%s.xtc'%self.npt.out]
            self.tprs = [self.npt.tpr]
        else:
            #self.sys_prepare(pdb)
            self.trajfiles = trajfiles
            self.tprs = tprs

    def sys_prepare(self, gro):

        # Generate topology file
        if self.top is None: self.top = 'topol.top'
        if ".gro" in gro:
            self.system.gen_top(gro=gro)
        else:
            self.system.gen_top()
        # Build box
        self.system.build_box()
        # Solvate the system
        self.system.solvate(top=self.top)
        # Add ions
        #if self.run['Cl']>0 or self.run['Na']>0:
        self.system.ionize(top=self.top)

    def sys_equilibrate(self, step, fn=None):

        if step is 'mdp':
            # Relax the system
            if fn is None: fn = 'min'
            self.min.gen_tpr(mdp=fn, tpr='data/minim.tpr')
            self.min.gmx_run(out='data/minim')
        elif step is 'nvt':
            # Obtain NVT ensemble
            if fn is None: fn = 'nvt'
            self.nvt.gen_tpr(mdp=fn, tpr='data/posre.tpr')
            self.nvt.gmx_run(out='data/posre')
        elif step is 'npt':
            # Obtain NPT ensemble
            if fn is None: fn = 'npt'
            self.npt.gen_tpr(mdp=fn, tpr='data/npt.tpr')
            self.npt.gmx_run(out='data/npt')

    def all_sys_equilibration(self):
        """ 
        Get at one time the whole thermodynamic ensemble

        """
        self.min = launcher.Launcher(self.system)
        self.sys_equilibrate('mdp')
        min_system = copy.deepcopy(self.min.system)
        self.nvt = launcher.Launcher(min_system)
        self.sys_equilibrate('nvt')
        nvt_system = copy.deepcopy(self.nvt.system)
        self.npt = launcher.Launcher(nvt_system)
        self.sys_equilibrate('npt')

    def adaptive_sampling(self, n_runs, lagt, mcs=85, ms=15, sym=False, rate_mat=True,    \
                        scoring='populations', n_epochs=2, nsteps=250000, max_time=100000.0,\
                        method=None, not_run=False):
        """
        Implementation of the Adaptive Sampling algorithm

        Parameters
        ----------
        n_runs : int
            Total number of parallel runs
        lagt : float
            Lag time (in nanosec)
        mcs : int
            min_cluster_size for HDBSCAN
        ms : int
            min_samples for HDBSCAN
        rate_mat : bool
            Compute K matrix, otherwise T will be calculated
        scoring : str
            Scoring function to use (str): count, popul
        n_epochs : int
            Number of outer loops in AS algorithm
        nsteps   : int
            Number of MD steps per epoch
        max_time : int
            Maximum simulation time in ps
        not_run : bool
            Use already simulated trajectories for epochs

        """
        n = 0
        offset = len(self.trajfiles)
        sim_time = 0.0
        while True:
            # ANALYZER #
            self.anal = analyzer.Analyzer(self.trajfiles)
            self.anal.build_msm(n, n_runs, lagt, method=method, mcs=mcs, ms=ms,\
                sym=sym, gro=self.gro_initial, rate_mat=rate_mat, offset=offset)
            if not_run:
                inputs = self.anal.resampler(self.tprs, scoring,\
                                        not_run=True)
                self.not_run(inputs,nsteps,n) #self.trajfiles
                n += 1
                if n > n_epochs: break
                continue
            else:
                inputs = self.anal.resampler(self.tprs, scoring)
            
            n += 1
            if n > n_epochs or sim_time > max_time:
#                # Call SuperMSM to converge last MSM
                break

            # PRODUCTION #
            tprs, outs = [], []
            i = 0
            for gro in inputs:
                i += 1
                self.system = system.System(gro, self.system.ff, self.system.wat,\
                            constraints=self.system.cons)#topology=self.top
                self.sys_prepare(gro)
                # Set up the simulation box and thermodynamic ensemble
                self.all_sys_equilibration()
                tpr, out = 'data/tpr/%s_%s.tpr'%(n, i), 'data/out/%s_%s'%(n, i)
                launcher_lib.checkfile(tpr)
                launcher_lib.checkfile(out)
                self.npt.gen_tpr(mdp='prod', tpr=tpr, nsteps=nsteps)#, sim_time=sim_time)
                tprs.append(tpr)
                outs.append(out)
                self.trajfiles.append('%s.xtc'%out)
                # CLEAN not interesting files!

            [self.tprs.append(tpr) for tpr in tprs]
            # Run parallel short trajectories
            self.npt.mp_run(tprs, outs)

    def not_run(self, inputs, nsteps, n):
        """
        Extract trajectory parts from h5 reference file according to 'inputs'

        """
        traj = self.anal.data[0]
        for inp in inputs:
            print(traj, n, inp[0], inp[1])
            data = traj[inp[0]:inp[0]+nsteps]
            # n stands for nepoch and inp[0] for label
            h5file = "./trajs/cossio_%g_%g.h5"%(n,inp[0])
            with h5py.File(h5file, "w") as trajectory:
                trajectory.create_dataset("data", data=data)
            self.trajfiles.append(h5file)

