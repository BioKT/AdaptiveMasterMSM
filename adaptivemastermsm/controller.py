"""
This file is part of the AdaptiveMasterMSM package.

"""
# Tools
import copy
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
    def __init__(self, pdb, forcefield=None, water=None,\
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
        self.sys_prepare(pdb)

        if trajfiles is None:
            self.all_sys_equilibration()
            self.trajfiles = ['%s.xtc'%self.npt.out]
            self.tprs = [self.npt.tpr]
        else:
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

    def adaptive_sampling(self, n_runs, lagt, mcs=185, ms=145, sym=False, rate_mat=True,\
                            scoring='populations', n_epochs=2, max_time=1000):
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
        max_time : int
            Maximum simulation time in ns

        """
        # Set options from 'equilibration' calculations
        forcefield = self.system.ff
        water = self.system.wat
        constraints = self.system.cons

        n = 0
        tprs = self.tprs
        while True:
            # ANALYZER #
            self.anal = analyzer.Analyzer(self.trajfiles)
            self.anal.build_msm(n, n_runs, lagt, method='ramagrid', \
                mcs=mcs, ms=ms, sym=sym, gro=self.gro_initial, rate_mat=rate_mat)
            inputs = self.anal.resampler(tprs, scoring=scoring)
            
            n += 1
            if n > n_epochs: # or sim_time > max_time:
#                # Call SuperMSM to converge last MSM
                break

            # PRODUCTION #
            tprs, outs, self.trajfiles = [], [], []
            i = 0
            for gro in inputs:
                i += 1
                self.system = system.System(gro, forcefield, water,\
                            constraints=constraints)#topology=self.top
                self.sys_prepare(gro)
                # Set up the simulation box and thermodynamic ensemble
                self.all_sys_equilibration()
                tpr, out = 'data/tpr/%s_%s.tpr'%(n, i), 'data/out/%s_%s'%(n, i)
                launcher_lib.checkfile(tpr)
                launcher_lib.checkfile(out)
                self.npt.gen_tpr(mdp='prod', tpr=tpr)
                tprs.append(tpr)
                outs.append(out)
                self.trajfiles.append('%s.xtc'%out)
                # CLEAN not interesting files!

            # Run parallel short trajectories
            self.npt.mp_run(tprs, outs)

