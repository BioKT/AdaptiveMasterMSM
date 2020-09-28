"""
This file is part of the AdaptiveMasterMSM package.

"""
# General
import h5py
import os, sys
from os import path
import itertools
# Maths
import numpy as np
# Plotting
import matplotlib.pyplot as plt
import seaborn as sns
sns.set(style="ticks", color_codes=True, font_scale=1.5)
sns.set_style({"xtick.direction": "in", "ytick.direction": "in"})
# Clustering
from sklearn.preprocessing import StandardScaler
import hdbscan
import mdtraj as md
# AdaptiveMasterMSM
from ..analyzer import analyzer_lib
# MasterMSM
from mastermsm.trajectory import traj
from mastermsm.msm import msm
from mastermsm.fewsm import fewsm

class Analyzer(object):
    """
    Build the MSM and resample to prepare next round
    """

    def __init__(self, trajfiles):
        """
        Create an updated MSM, resample, and generate new input

        Parameters
        ----------
        trajfiles : str
            Path to trajectory files separated by spaces only

        """
        # Read parallel trajectories
        self.data = []
        print(trajfiles.split())
        for f in trajfiles.split():
            if not path.isfile(f):
                raise ValueError("Trajectory file(s) not valid")
            if '.xtc' in f:
                self.data.append(f)
            elif '.h5' in f:
                fr = h5py.File(f, 'r')
                self.data.append(np.array(fr['data']))
                fr.close()

        # Next steps better to do from outside:
        # Do MSM and get macrostates build_msm(mcs, ms, lagt, rate_mat=True)
        # Call resampler(macmsm) and gen_input

        # Note: gen_input is only called from outside (controller)

    def build_msm(self, n_runs, lagt, dt=1, mcs=185, ms=145, sym=True,\
                    n_clusters=0, rate_mat=True, gro=None):
        """
        Build the MSM for the next round.
        Return macrostates obtained from the MSM.
        
        Parameters
        ----------
        n_runs : int
            Number of parallel simulations
        lagt : float
            Lag time (in nanosec)
        mcs : int
            min_cluster_size for HDBSCAN
        ms : int
            min_samples for HDBSCAN
        dt : float
            Time step (in ns) between subsequent frames in trajectory
        sym : bool
            Impose symmetry on Count matrix, i.e. impose detailed balance
        n_clusters : int
            Number of clusters or macrostates
        rate_mat : bool
            Build the MSM from K matrix, otherwise T is used
        gro : str
            Path to .gro or .pdb file

        """
        self.n_runs, self.n_clusters = n_runs, n_clusters
        self.lt, self.dt, self.sym, self.rate = lagt, dt, sym, rate_mat
        self.min_cluster_size, self.min_samples = mcs, ms

        #1- Clustering:
        # one 'labels' per parallel run, TimeSeries will merge all together
        self.labels_all = []
        trs = []
        for i in range(len(self.data)):
            #labels = self.gen_clusters_mueller(i)
            labels, tr = self.gen_clusters_ramachandran(i, gro, method='hdbscan',\
                            mcs=self.min_cluster_size, ms=self.min_samples)
            trs.append(tr)
            self.labels_all.append(labels)
        
        self.trajs = trs

        #2- do MSM:
        #self.gen_msm(tr_instance=False)
        self.gen_msm(tr_instance=True)

        if self.n_clusters > 0:
            self.gen_mmsm()
        else:
            self.mMSM = self.MSM
            self.n_clusters = len(self.mMSM.keep_states)

        #3- Wrap up
        print("MSM done")
        sys.stdout.flush()

    def gen_clusters_ramachandran(self, i, gro, method='ramagrid', mcs=185, ms=185):
        """
        Cluster trajectories into microstates by using Ramachandran angle regions
        Return labels (int array) containing trajectory by microstates

        """
        #phi = md.compute_phi(tr.mdt)
        #psi = md.compute_psi(tr.mdt)
        #tr = traj.TimeSeries(top=gro, traj=[self.data[i]])
        tr = traj.TimeSeries(top=gro, traj=self.data[i])
        if method == 'ramagrid':
            tr.discretize(method='ramagrid', nbins=5)
        elif method == 'hdbscan':
            mcs, ms = self.min_cluster_size, self.min_samples
            tr.discretize(method='hdbscan', mcs=mcs, ms=ms)
        else:
            tr.discretize(states=['A', 'E'])
        
        tr.find_keys()
        tr.keys.sort()
        """
        fig, ax = plt.subplots(figsize=(10,3))
        ax.plot(tr.mdt.time, [x for x in tr.distraj], lw=2)
        ax.set_xlim(0,5000)
        ax.set_ylim(0.8,2.2)
        ax.set_xlabel('Time (ps)', fontsize=20)
        ax.set_ylabel('state', fontsize=20)
        plt.savefig('traj.png')
        """
        return tr.distraj, tr

    def gen_clusters_mueller(self, i):
        """
        Cluster trajectories into microstates by using HDBSCAN.
        Return labels (int array) containing trajectory by microstates

        """
        data = self.data[i]
        X = data[:,[1,2]]
        hb = hdbscan.HDBSCAN(min_cluster_size = self.min_cluster_size, \
                            min_samples = self.min_samples).fit(X)
        labels = hb.labels_
        n_micro_clusters = len(set(labels)) - (1 if -1 in labels else 0)
        n_noise = list(labels).count(-1)

        # remove from clusters points with small (<0.1) probability
        for i, x_i in enumerate(labels):
            if hb.probabilities_[i] < 0.1:
                labels[i] = -1

        # # plot better from outside, put plotting function in analyzer_lib
        # colors = ['royalblue', 'maroon', 'forestgreen', 'mediumorchid', \
        # 'tan', 'deeppink', 'olive', 'goldenrod', 'lightcyan', 'lightgray']
        # vectorizer = np.vectorize(lambda x: colors[x % len(colors)])
        # plt.scatter(X[:,0],X[:,1], c=vectorizer(labels))
        # plt.savefig('mueller_hdbscan.png')

        # remove noise from microstate trajectory
        i = 0
        last = labels[i]
        while last == -1:
            i += 1
            last = labels[i]

        for i, x_i in enumerate(labels):
            if x_i == -1:
                labels[i] = last
            else:
                last = x_i

        return labels

    def gen_msm(self, tr_instance=False):
        """
        Do MSM and build macrostates

        Parameters
        ----------
        tr_instance : bool
            True if an instance of trajectory exists

        """

        # Create an instance of traj for each trajectory and merge in a single list
        if tr_instance:
            distrajs = self.trajs
        else:
            labels = self.labels_all
            dt = self.dt
            distrajs = []
            for label in labels:
                distraj = traj.TimeSeries(distraj=list(label), dt=dt)
                distraj.find_keys()
                distraj.keys.sort()
                distrajs.append(distraj)

        lagt = self.lt
        sym = self.sym
        # Invoke SuperMSM to collect keys from all trajectories
        smsm = msm.SuperMSM(distrajs, sym=sym)
        micro_msm = msm.MSM(data=distrajs, keys=smsm.keys, lagt=lagt, sym=sym)
        micro_msm.do_count()
        if self.rate:
            #smsm = msm.SuperMSM(distrajs, sym=True) #keys=distraj.keys
            smsm.do_lbrate()
            micro_msm.do_rate(method='MLPB', evecs=False, init=smsm.lbrate)
        else:
            micro_msm.do_trans(evecs=False)
        
        fig, ax = plt.subplots()
        #ax.errorbar(range(1,2),np.log(micro_msm.tauK[0:2]), fmt='o-')
        ax.errorbar(range(1,2),micro_msm.tauT[0:10], fmt='o-')
        ax.set_xlabel('Eigenvalue index')
        ax.set_ylabel(r'$\tau_i$ (ns)')
        plt.savefig('eigs_T.png')
        
        self.MSM = micro_msm

    def gen_mmsm(self):
        """
        Build macrostates

        """

        n_clusters = self.n_clusters
        MSM = self.MSM
        mmsm = fewsm.FEWSM(MSM, N=n_clusters)
        self.mMSM = mmsm

        import matplotlib.cm as cm
        fig, ax = plt.subplots(figsize=(5,5))
        mat = np.zeros((20,20), float)
        for i in MSM.keep_keys:
            j = MSM.keep_keys.index(i)
            if j in mmsm.macros[0]:
                mat[i%20, int(i/20)] = 1
            elif j in mmsm.macros[1]:
                mat[i%20, int(i/20)] = 2
            else:
                mat[i%20, int(i/20)] = 3
            mat#print i, i[0]%20, int(i[0]/20), -i[1]
        my_cmap = cm.get_cmap('viridis')
        my_cmap.set_under('w')
        ax.imshow(mat.transpose(), interpolation="none", origin='lower', \
             cmap=my_cmap, vmin = 0.5)
        plt.savefig('fewms.png')
        
    def resampler(self, scoring='populations'):
        """
        
        Parameters
        ----------
        scoring : str
            Scoring function to determine how to resample

        """

        if scoring == "counts":
            states = self.counts()
        elif scoring == "populations":
            states = self.populs()
        elif scoring == "non_detailed_balance":
            if self.sym: raise Exception("Cannot impose symmetry with chosen scoring criteria.")
            states = self.non_detbal()

        self.gen_input(states)

    def non_detbal(self):
        """
        Resampling criteria proportional to lack of detailed
        balance of each MSM state

        """

        nondb = abs( self.MSM.count - np.transpose(self.MSM.count) )
        return np.sum(nondb, axis=1)        

    def populs(self):
        """
        Resampling criteria inversely proportional to
        MSM's populations

        """

        return 1./self.MSM.peqK if self.rate else 1./self.MSM.peqT

    def counts(self):
        """
        Resampling criteria inversely proportional to
        how many times those states have been visited

        """

        return 1./np.sum(self.MSM.count, axis=1)

    def gen_input(self, states):
        """
        
        Parameters
        ----------
        inpcrd : str
            File containing coordinates for next round
        states : np array
            Scoring of each macrostate
        n_msm_runs : int array
            Number of runs for each macrostate

        """
        # Determine distribution of new runs according to 'states'
        n_runs = self.n_runs
        n_msm_runs = []
        for x_i in states:
            if np.log(x_i) > 0.5:
                n_msm_runs.append(np.log(x_i))
            else:
                n_msm_runs.append(0.0)
        fac = n_runs / np.sum(n_msm_runs)
        n_msm_runs = [fac*x_i for x_i in n_msm_runs]
        n_msm_runs = [round(x_i) for x_i in n_msm_runs]
        
        print(len(n_msm_runs), np.sum(n_msm_runs))
        print(n_msm_runs, n_runs)
        
        # Use weights to randomly choose a frame and create a corresponding .gro file
        analyzer_lib.gen_input_weights(n_msm_runs, self.labels_all, self.trajs)

        # Use a dict containing all trajs, labels and frames to pick randomly a frame
        # to do ...


