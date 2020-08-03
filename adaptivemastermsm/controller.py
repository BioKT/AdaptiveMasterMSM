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
from adaptivemastermsm.launcher import launcher
###from adaptivemastermsm.analyzer import analyzer
###from adaptivemastermsm.system import system

class Controller(object):
    """
    Driver of the adaptive sampling algorithm
    """

    def __init__(self, pdb_fn, forcefield, water, nproc=1):
        """
        Iteratively, create a GROMACS file, launch, and analyze.

        Args:
            MD:
                forcefield (str): amber96, ...
                water: water model for GROMACS (tip3p, ...)
            MSM:
                mcs (int): min_cluster_size for HDBSCAN
                ms (int): min_samples for HDBSCAN
                lagt (float): lag time (in nanosec)
                rate (bool): Compute K matrix, otherwise T will be calculated
            AS:
                pdb_fn (str): .pdb filename, also .gro can be provided
                scoring function to use (str): count, popul
                n_rounds (int): number of outer loops in AS algorithm
                max_time (int): maximum simulation time in ns
                n_samp_clus(int): number of samples per cluster
        """

        # GLOBAL PARAMETERS #
        self.pdb_fn = pdb_fn
        self.ff = forcefield
        self.water = water

        n_rounds = 10
        max_time = 100
        
        # EQUILIBRATION #
        lau_eq = launcher.Launcher('Equilibration', nproc, self.ff, self.water, \
        self.pdb_fn, 'wetfn', dry_xtc_file='dryfn', last_wet_snapshot='lastwetfn')

        n = 0
        while True:
            n += 1
            # somewhere here do MSM
            # analyzer.Analyzer()

            # PRODUCTION #
            # launcher.Launcher()

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

            if (n > n_rounds or sim_time > max_time):
                # Call SuperMSM to converge last MSM
                break
        
