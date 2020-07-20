"""
This file is part of the MasterMSM package.

"""
# Tools
import h5py
import random
rand = lambda: random.randint(0, 255)
# Maths
import numpy as np
# Plotting
import matplotlib.pyplot as plt
# Clustering
from sklearn.preprocessing import StandardScaler
import hdbscan
# MasterMSM
from mastermsm.trajectory import traj
from mastermsm.msm import msm

class run_msm(object):
    """
    Build the MSM from short trajectories and update new clusters
    """

    def __init__(self, trajfile, **kwards):
        """
        From a new round, create an updated MSM

        Args:
            trajfile (str): File to read trajectories
        """

        # Read trajectory from controller preprocessor
        self.h5file = trajfile
        f = h5py.File(self.h5file, 'r')
        self.data = np.array(f['data'])
        f.close()


    def builder(self):
        """
        Build the MSM for the next round
        """
        #1- clustering: gen_clusters

        #2- MSM(data, lagt, sym)--> do_count()--> do_rate or do_trans
        distraj = traj.TimeSeries(distraj=list(self.labels), dt=1)
        distraj.find_keys()
        distraj.keys.sort()

        Could we skip SuperMSM?


        #3- if resample: save clusters into file and call to input_gen()
        #   else: finish

    def gen_clusters(self, msc=50, ms=50, plotter):
        """
        Cluster trajectories into microstates

        Args:
            data: np array containing trajectory
            msc (int): min_cluster_size for HDBSCAN
            ms (int): min_samples for HDBSCAN
            plotter (boolean): plot microstates or not
        """
        X = data[:,[1,2]] #data[:1000,[1,2]]
        hb = hdbscan.HDBSCAN(min_cluster_size=185,min_samples=145).fit(X)
        self.labels = hb.labels_
        n_clusters_ = len(set(self.labels)) - (1 if -1 in self.labels else 0)
        n_noise_ = list(self.labels).count(-1)

        # exclude from clusters points with small (<0.1) probability
        for i, x_i in enumerate(self.labels):
            if hb.probabilities_[i] < 0.1:
                self.labels[i] = -1

        # remove noise from microstate trajectory
        last = self.labels[0]
        for i, x_i in enumerate(self.labels):
            if x_i == -1:
                self.labels[i] = last
            else:
                last = x_i
        self.labels[0] = self.labels[1]

        call to plotter if requested (tree_plot and/or self.labels)

