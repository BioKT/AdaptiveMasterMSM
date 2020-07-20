"""
This file is part of the AdaptiveMasterMSM package.

"""
# Tools
import h5py
import random
rand = lambda: random.randint(0, 255)
# Maths
import numpy as np
# Plotting
import matplotlib.pyplot as plt

class controller(object):
    """
    Driver of the adaptive sampling algorithm
    """

    def __init__(self, asfile):
        """
        Iteratively, create a GROMACS file, launch, and analyze.

        Args:
            asfile (str): file to read general information
            MD: molecule, model (GROMACS)
            MSM:
                mcs (int): min_cluster_size for HDBSCAN
                ms (int): min_samples for HDBSCAN
                lagt (float): lag time (in nanosec)
                rate (bool): Compute K matrix, otherwise T will be calculated
            AS:
                scoring function to use: count, popul
                n_rounds (int): number of outer loops in AS algorithm
                n_samp_clus(int): number of samples per cluster
        """

        self.asfile = asfile
        self.asfile = ConfigParser(interpolation=None)
        self.asfile.read(self.asfile)
        # solo para que vea un ejemplo (mirar sampler.example.cfg)
        self.featurized = self.config.get("system", "featurized",
                                          fallback="featurized_%d")
        # Input general
        self.epoch = int(kwargs.get("generation",
                        self.config.getint("production", "generation",
                                           fallback=0)))
        # Input for HDBSCAN
        self.mcs = 
        self.ms  = 
        # Input for MSM
        self.lagt = 
        self.rate = 

    def driver(self):
        """
        Calls 3 classes: system, launcher, and analyzer
        """

        #1- system

        #2- launcher

        #3- analyzer(trajfile, self.epoch, self.mcs, self.ms, \
        #           self.lagt, self.rate)

    def shell_out(self, commands):
        """ executes commands intelligently via subprocess
            - waits for serial completion 
        """
    
        log = open('driver.log','w')
        for command in commands:
            print "Executing: %s" % command
            x = subprocess.Popen( command, bufsize=-1, shell=True, stdout=log, stderr=log )
            x.wait()
            if x.returncode != 0:
                raise Exception("Error: '%s' quit with exit status: %d" % (command, x.returncode) )
        log.close()
    
        return

    def clean_working_directory(self):

        print "Cleaning up..."
        os.system('rm \#*')
        os.system('mkdir logs')
        os.system('mv *.log logs')
        os.system('rm *.trr *.top *.itp *.edr *.cpt *.tpr out.gro conf.gro')
        os.system('mv equilibration.xtc *.mdp *.gro logs')
    
        return

