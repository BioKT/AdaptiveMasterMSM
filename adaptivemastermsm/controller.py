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

        """
        # try for instance pdb_fn='alaTB.gro', forcefield='amber96', water='tip3p'
        self.pdb_fn = pdb_fn
        self.ff = forcefield
        self.water = water        
        # EQUILIBRATION #
        lau_eq = launcher.Launcher('Production', self.ff, self.water, \
                self.pdb_fn, 'wetfn', dry_xtc_file='dryfn', \
                last_wet_snapshot='lastwetfn')
        """


    def runner(self, md_step):

        self.md_step = md_step
        # Take production or equilibration paths
        if self.md_step == 'Equilibration':
            self.driver_equilibration()
        elif self.md_step == 'Production':
            self.driver_production()
        else:
            print ("md_step %s not valid" % self.md_step)
            raise Exception("It must be 'Production' or 'Equilibration'")

        return
        
    def driver_equilibration(self):
    
        equil = self.set_equilibration()
        self.run = {}
        equil['total_steps'] = \
            int(equil['total_time'] / equil['timestep'])
        
        for p in equil.keys():
            self.run[p] = equil[p]

        return           

    def driver_production(self):

        prod = self.set_production()
        self.run = {}
        prod['total_steps'] = \
            int(prod['total_time'] / prod['timestep'])

        for p in prod.keys():
            self.run[p] = prod[p]

        return

    def set_production(self):
        """
        Parameters for the MD production run

        """
        production = {
            'timestep'       : 0.002,         # ps
            'total_time'     : 1,             # ps / 10ns
            'log_freq'       : 10000,         # timesteps / 20ps
            'xtc_freq'       : 100,           # timesteps
            'temperature'    : 300,           # Kelvin
            'thermostat'     : 'Nose-Hoover', # type
            'box_size'       : 3.5,           # Angstroms
            'barostat'       : 'No',          # No barostat
            'pressure'       : 1,             # pressure
            'Cl'             : 0,             # number Cl's
            'Na'             : 0              # number Na's
        }

        return production

    def set_equilibration(self):
        """
        Parameters for the MD equilibration run

        """
        equilibration = {
            'timestep'       : 0.002,         # ps
            'total_time'     : 1,             # ps
            'barostat'       : 'berendsen',   # type
            'pressure'       : 1,             # bar
            'log_freq'       : 10000,         # timesteps / 20ps
            'xtc_freq'       : 100,           # timesteps
            'temperature'    : 300,           # Kelvin
            'thermostat'     : 'Nose-Hoover', # type
            'box_size'       : 1.0,           # Angstroms
            'Cl'             : 0,             # number Cl's
            'Na'             : 0              # number Na's
        }

        return equilibration

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

            
        

