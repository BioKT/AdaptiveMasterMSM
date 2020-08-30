"""
This file is part of the AdaptiveMasterMSM package.

"""
# Tools
import h5py
###import random
###rand = lambda: random.randint(0, 255)
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
    def __init__(self):
        """
        Iteratively, create a GROMACS file, launch, and analyze.

        Parameters
        ----------
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
        # GLOBAL PARAMETERS #
        

    def run_parameters(self, md_step):

        self.md_step = md_step
        self.run = {}
        # Take production or equilibration paths
        if self.md_step == 'Equilibration':
            controller_lib.driver_equilibration(self.run)
        elif self.md_step == 'Production':
            controller_lib.driver_production(self.run)
        else:
            print ("md_step %s not valid" % self.md_step)
            raise Exception("It must be 'Production' or 'Equilibration'")

        return

    def sys_equilibration(self, gro, forcefield=None, water=None,\
        topology=None, constraints=None):

        pdb = system.System(gro, forcefield=forcefield, water=water,\
                 topology=topology, constraints=constraints)
        # Generate topology file
        if topology is None:
            topology='topol.top'
            if ".gro" in gro:
                pdb.gen_top(gro=gro, top=topology)
            else:
                pdb.gen_top(top=topology)
        # Build box
        pdb.build_box(gro=gro, d=self.run['box_size'])
        # Solvate the system
        pdb.solvate(top=topology)
        # Add ions
        inp = 'conf_solv.gro'
        if self.run['Cl']>0 or self.run['Na']>0:
            pdb.ionize(anion=self.run['Cl'], cation=self.run['Na'],\
                        top=topology)
            inp = 'conf_solv_ions.gro'
        # Equilibrate the system, and obtain stable NVT-NPT ensemble
        pdb.minimize(inp=inp, top=topology)
        pdb.nvt(mdp='nvt.mdp', top=topology)
        pdb.npt(mdp='npt.mdp', top=topology)


    def sys_production(self, gro, forcefield, water,\
                        gro_files=None, top_files=None):

        pdb = launcher.Launcher(self.md_step, gro, forcefield, water)
        # Generate input files for MDRUN
        if gro_files is None: gro_files='npt.gro'
        if top_files is None: top_files='topol.top'
        gmxinp = pdb.inputs(gro_files, top_all=top_files,\
                    chk='npt.cpt', mdp_params=self.run)
        # Run simulation using Gromacs
        pdb.runner(gmxinp)

    def adaptive_sampling(self):
        
        n_rounds = 10
        max_time = 100
        n = 0
        while True:
            n += 1
            # somewhere here do MSM
            # analyzer.Analyzer()

            # PRODUCTION #
            # launcher.Launcher()

            if n > n_rounds or sim_time > max_time:
                # Call SuperMSM to converge last MSM
                break

        """
        # define multiprocessing options
        nproc = mp.cpu_count()
        pool = mp.Pool(processes=nproc)
        # generate multiprocessing input
        sliding = True
        mpinput = [[x.distraj, x.dt, self.keys, x.dt, sliding]
                for x in self.data]
        # run counting using multiprocessing
        result = pool.map(msm_lib.calc_count_worker, mpinput)
        pool.close()
        pool.join()
        # add up all independent counts
        count = reduce(lambda x, y: np.matrix(x) + np.matrix(y), result)
        """

            
        

