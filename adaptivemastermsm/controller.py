"""
This file is part of the AdaptiveMasterMSM package.

"""
# Tools
import copy
#import h5py
#import random
#rand = lambda: random.randint(0, 255)
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

class Controller(object):
    """
    Driver of the adaptive sampling algorithm

    """
    def __init__(self, pdb, forcefield=None, water=None,\
                topology=None, constraints=None):
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
        mcs : int
            min_cluster_size for HDBSCAN
        ms : int
            min_samples for HDBSCAN
        lagt : float
            lag time (in nanosec)  
        rate : bool
            Compute K matrix, otherwise T will be calculated
        scoring : str
            Scoring function to use (str): count, popul
        n_rounds : int
            Number of outer loops in AS algorithm
        max_time : int
            Maximum simulation time in ns
        n_samp_clus : int
            Number of samples per cluster

        """

        # Define the system
        Controller.system = system.System(pdb, forcefield, water,\
                            topology=topology, constraints=constraints)
        self.top = topology
        self.sys_prepare(pdb)
        
        # Prepare the simulation box and thermodynamic ensemble
        Controller.min = launcher.Launcher(Controller.system)
        self.sys_equilibrate('mdp')
        min_system = copy.deepcopy(Controller.min.system)
        Controller.nvt = launcher.Launcher(min_system)
        self.sys_equilibrate('nvt')
        nvt_system = copy.deepcopy(Controller.nvt.system)
        Controller.npt = launcher.Launcher(nvt_system)
        self.sys_equilibrate('npt')

        # Production
        npt_system = copy.deepcopy(Controller.npt.system)

    def sys_prepare(self, gro):

        # Generate topology file
        if self.top is None:
            topology = 'topol.top'
            if ".gro" in gro:
                self.system.gen_top(gro=gro)
            else:
                self.system.gen_top()
        # Build box
        self.system.build_box()
        # Solvate the system
        self.system.solvate(top=topology)
        # Add ions
        #if self.run['Cl']>0 or self.run['Na']>0:
        self.system.ionize(top=topology)

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

#    def sys_production(self, gro, forcefield, water,\
#                        gro_files=None, top_files=None):
#
#        pdb = launcher.Launcher(self.md_step, gro, forcefield, water)
#        # Generate input files for MDRUN
#        if gro_files is None: gro_files='npt.gro'
#        if top_files is None: top_files='topol.top'
#        gmxinp = pdb.inputs(gro_files, top_all=top_files,\
#                    chk='npt.cpt', mdp_params=self.run)
#        # Run simulation using Gromacs
#        pdb.runner(gmxinp)
#
#    def run_parameters(self, md_step):
#
#        self.md_step = md_step
#        self.run = {}
#        # Take production or equilibration paths
#        if self.md_step == 'Equilibration':
#            controller_lib.driver_equilibration(self.run)
#        elif self.md_step == 'Production':
#            controller_lib.driver_production(self.run)
#        else:
#            print ("md_step %s not valid" % self.md_step)
#            raise Exception("It must be 'Production' or 'Equilibration'")
#
#    def adaptive_sampling(self):
#        
#        n_rounds = 10
#        max_time = 100
#        n = 0
#        while True:
#            n += 1
#            # somewhere here do MSM
#            # analyzer.Analyzer()
#
#            # PRODUCTION #
#            # launcher.Launcher()
#
#            if n > n_rounds or sim_time > max_time:
#                # Call SuperMSM to converge last MSM
#                break
#
#        """
#        # define multiprocessing options
#        nproc = mp.cpu_count()
#        pool = mp.Pool(processes=nproc)
#        # generate multiprocessing input
#        sliding = True
#        mpinput = [[x.distraj, x.dt, self.keys, x.dt, sliding]
#                for x in self.data]
#        # run counting using multiprocessing
#        result = pool.map(msm_lib.calc_count_worker, mpinput)
#        pool.close()
#        pool.join()
#        # add up all independent counts
#        count = reduce(lambda x, y: np.matrix(x) + np.matrix(y), result)
#        """

            
        

