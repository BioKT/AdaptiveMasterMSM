"""
This file is part of the AdaptiveMasterMSM package.

"""
# General
import h5py
import sys
# Maths
import numpy as np
# Clustering
from sklearn.preprocessing import StandardScaler
import hdbscan
# MasterMSM
from mastermsm.trajectory import traj
from mastermsm.msm import msm

class analyzer(object):
    """
    Build the MSM and resample to prepare next round
    """

    def __init__(self, trajfile, epoch, mcs, ms, lagt, rate):
        """
        Create an updated MSM, resample, and generate new input

        Args:
            trajfile (str): File to read trajectories
            epoch (int): Which epoch/round/generation we are at
            mcs (int): min_cluster_size for HDBSCAN
            ms (int): min_samples for HDBSCAN
            lagt (float): lag time (in nanosec)
            rate (bool): Compute K matrix, otherwise T will be calculated
        """

        # Read trajectory from controller preprocessor
        h5file = trajfile
        f = h5py.File(h5file, 'r')
        self.data = np.array(f['data'])
        f.close()

        # self.generation = epoch
        # self.min_cluster_size, self.min_samples = mcs, ms
        # self.lt, self.rate = lagt, rate

        # Call macmsm = build_msm()

        # Call resampler(macmsm) if self.generation < some maximum

        # Note: gen_input is only called from outside (controller)

    def build_msm(self):
        """
        Build the MSM for the next round
        
        Args:
            data: np array containing trajectories

        Return:
            macmsm: macrostates obtained from the MSM
        """
        print("EPOCH %d:" % self.generation)

        #1- clustering: labels = gen_clusters()

        #2- MSM: gen_msm(labels)

        #3- if self.generation > some maximum: Call SuperMSM to converge last MSM
        #   else: save clusters into file and finish this routine

        print("MSM done")
        sys.stdout.flush()

        return macmsm

    def gen_clusters(self):
        """
        Cluster trajectories into microstates by using HDBSCAN

        Return:
            labels (int array): trajectory by microstates
        """
        X = self.data[:,[1,2]]
        hb = hdbscan.HDBSCAN(min_cluster_size = self.min_cluster_size, \
                            min_samples = self.min_samples).fit(X)
        labels = hb.labels_
        n_clusters_ = len(set(labels)) - (1 if -1 in labels else 0)
        n_noise_ = list(labels).count(-1)

        # remove from clusters the points with small (<0.1) probability
        for i, x_i in enumerate(labels):
            if hb.probabilities_[i] < 0.1:
                labels[i] = -1

        # remove noise from microstate trajectory
        last = labels[0]
        for i, x_i in enumerate(labels):
            if x_i == -1:
                labels[i] = last
            else:
                last = x_i
        labels[0] = labels[1]

        # plot if requested (tree_plot and/or labels)

        return labels


    def gen_msm(self, labels):
        """
        Do MSM and build macrostates

        Args:
            labels (int array): trajectory by clusters
        """

        #1- msm = MSM(data=labels, lagt=self.lt, sym=False)

        #2- msm.do_count()

        #3- if self.rate: msm.do_rate(method='Taylor', evecs=False, init=False)
        #   else: msm.do_trans(evecs=False)

        #4- PCCA+?
        
    def resampler(self, macmsm):
        """
        Args:
            macmsm: macrostates obtained from the MSM

        Return:
            Information containing how to resample
        """

    def gen_input(self):
        """
        Return:
            inpcrd (str): File containing coordinates for
                next round, read by class system
        """
