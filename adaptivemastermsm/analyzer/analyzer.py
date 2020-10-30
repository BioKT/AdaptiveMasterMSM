"""
This file is part of the AdaptiveMasterMSM package.

"""
# General
import h5py
import os, sys
from os import path
import random
# Maths
import numpy as np
# Plotting
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import seaborn as sns
sns.set(style="ticks", color_codes=True, font_scale=1.5)
sns.set_style({"xtick.direction": "in", "ytick.direction": "in"})
# Clustering
#from sklearn.preprocessing import StandardScaler
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
        trajfiles : list
            Path(s) to trajectory file(s)

        """
        # Read parallel trajectories
        self.data = []
        for f in trajfiles:
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

    def build_msm(self, n_epoch, n_runs, lagt, mcs=185, ms=145, dt=1, sym=True,\
                    n_clusters=0, rate_mat=True, gro=None, method='hdbscan'):
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
        method : str
            Discretization method (hdbscan, rama, ramagrid)

        """
        self.n_epoch, self.n_runs, self.n_clusters = n_epoch, n_runs, n_clusters
        self.lt, self.dt, self.sym, self.rate = lagt, dt, sym, rate_mat
        self.min_cluster_size, self.min_samples = mcs, ms
        #1- Clustering:
        # one 'labels' per parallel run, TimeSeries will merge all together
        self.labels_all = []
        trs = []
        for i in range(len(self.data)):
            #labels = self.gen_clusters_mueller(i)
            labels, tr = self.gen_clusters_ramachandran(i, gro, method,\
                            self.min_cluster_size, self.min_samples)
        
            data = np.column_stack((tr.mdt.time, [x for x in tr.distraj]))
            h5file = "data/out/%g_%g_traj.h5"%(self.n_epoch, i)
            with h5py.File(h5file, "w") as hf:
                hf.create_dataset("rama_trajectory", data=data)

            trs.append(tr)
            self.labels_all.append(labels)
        self.trajs = trs

        #2- do MSM:
        #self.gen_msm(tr_instance=False)
        self.gen_msm(tr_instance=True)

        if self.n_clusters > 0:
            self.gen_macro_msm()
        else:
            self.mMSM = self.MSM
            self.n_clusters = len(self.mMSM.peqK) if self.rate else len(self.mMSM.peqT)

        #3- Wrap up
        print("MSM done")
        sys.stdout.flush()

    def gen_clusters_ramachandran(self, i, gro, method, mcs, ms):
        """
        Cluster trajectories into microstates by using Ramachandran angle regions
        Return labels (int array) containing trajectory by microstates

        """
        tr = traj.TimeSeries(top=gro, traj=self.data[i])#traj=[self.data[i]]
        if method == 'ramagrid':
            phi, psi = tr.discretize(method='ramagrid', nbins=20)
        elif method == 'hdbscan':
            mcs, ms = self.min_cluster_size, self.min_samples
            phi, psi = tr.discretize(method='hdbscan', mcs=mcs, ms=ms)
        else:
            phi, psi = tr.discretize(states=['A', 'E'])
        
        tr.find_keys()
        tr.keys.sort()

        data = np.column_stack((phi[1], psi[1], tr.distraj))
        h5file = "data/out/%g_%g_traj_ramachandran.h5"%(self.n_epoch, i)
        with h5py.File(h5file, "w") as hf:
            hf.create_dataset("rama_angles", data=data)

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
        lagt, dt, sym, rate = self.lt, self.dt, self.sym, self.rate
        # Create an instance of traj for each trajectory and merge in a single list
        if tr_instance:
            distrajs = self.trajs
        else:
            labels = self.labels_all
            distrajs = []
            for label in labels:
                distraj = traj.TimeSeries(distraj=list(label), dt=dt)
                distraj.find_keys()
                distraj.keys.sort()
                distrajs.append(distraj)
        # Invoke SuperMSM to collect keys from all trajectories
        smsm = msm.SuperMSM(distrajs, sym=sym)
        micro_msm = msm.MSM(data=distrajs, keys=smsm.keys, lagt=lagt, sym=sym)
        micro_msm.do_count()
        if rate:
            #smsm = msm.SuperMSM(distrajs, sym=True) #keys=distraj.keys
            smsm.do_lbrate()
            micro_msm.do_rate(method='MLPB', evecs=False, init=smsm.lbrate)
        else:
            micro_msm.do_trans(evecs=True)
        self.MSM = micro_msm

#    def gen_macro_msm(self, n_clusters=1):
#        """
#        Build macrostates
#
#        """
#        MSM = self.MSM
#        macro_msm = fewsm.FEWSM(MSM, N=n_clusters)
#        self.mMSM = macro_msm
#
#        fig, ax = plt.subplots(figsize=(5,5))
#        mat = np.zeros((20,20), float)
#        for i in MSM.keep_keys:
#            j = MSM.keep_keys.index(i)
#            if j in mmsm.macros[0]:
#                mat[i%20, int(i/20)] = 1
#            elif j in mmsm.macros[1]:
#                mat[i%20, int(i/20)] = 2
#            else:
#                mat[i%20, int(i/20)] = 3
#            mat#print i, i[0]%20, int(i[0]/20), -i[1]
#        my_cmap = cm.get_cmap('viridis')
#        my_cmap.set_under('w')
#        ax.imshow(mat.transpose(), interpolation="none", origin='lower', \
#             cmap=my_cmap, vmin = 0.5)
#        plt.savefig('fewms.png')
        
    def resampler(self, tprs, scoring='populations'):
        """
        
        Parameters
        ----------
        tprs : list
            tpr files from Launcher to generate new gro files
        scoring : str
            Scoring function to determine how to resample

        """

        if scoring == "counts":
            states = self.counts()
            inputs = self.gen_input(states.real, tprs)
        elif scoring == "populations":
            states = self.populs()
            inputs = self.gen_input(states.real, tprs)
        elif scoring == "non_detailed_balance":
            if self.sym:
                raise Exception("Cannot impose symmetry with chosen scoring criteria.")
            states = self.non_detbal()
            inputs = self.gen_input(states, tprs)
        elif scoring == "flux":
            states = self.flux_inbalance()
            inputs = self.gen_input(states, tprs)

        return inputs

    def populs(self):
        """
        Resampling probability inversely proportional to
        state populations

        """
        p = 1./self.MSM.peqT
        return p/np.sum(p)

    def counts(self):
        """
        Resampling probability inversely proportional to
        number of visits

        """
        #p = 1./np.sum(self.MSM.count, axis=1)
        p = [(1./np.sum(self.MSM.count[x, :])) \
            if np.sum(self.MSM.count) != 0     \
            else 0                             \
            for x in self.MSM.keep_states]
        if np.sum(p) == 0:
            sys.exit(" Error in 'resampler'. Please choose another 'scoring' option")
        return p/np.sum(p)

    def non_detbal(self):
        """
        Resampling probability proportional to lack of detailed
        balance for each MSM state

        """
        total = self.MSM.count + np.transpose(self.MSM.count)
        nondb = abs(self.MSM.count - np.transpose(self.MSM.count))/total
        nondb = np.sum(np.nan_to_num(nondb, 0), axis=1)
        # remove all states not satisfying ergodicity
        states = [nondb[i] for i in self.MSM.keep_states]
        states /= np.sum(states)
        print(len(self.MSM.keep_states), len(states))
        return states

    def flux_inbalance(self):
        """
        Resampling probality proportional to flux imbalance

        """
        #flux = [abs(np.sum(self.MSM.count[x, :]) - np.sum(self.MSM.count[:,x]))/\
        #        (np.sum(self.MSM.count[:,x]) + np.sum(self.MSM.count[x, :])) \
        #        for x in self.MSM.keep_states]
        flux = [(abs(np.sum(self.MSM.count[x,:]) - np.sum(self.MSM.count[:,x]))/\
        (np.sum(self.MSM.count[:,x]) + np.sum(self.MSM.count[x,:]))) \
        if (np.sum(self.MSM.count[:,x]) + np.sum(self.MSM.count[x,:])) >= 1
        else 0 \
        for x in self.MSM.keep_states] #keep_keys
        #this is done by above abs: flux += np.min(flux)
        if np.sum(flux) == 0:
            sys.exit(" Error in 'resampler'. Please choose another 'scoring' option")
        return flux/np.sum(flux)

    def gen_input(self, states, tprs):
        """
        
        Parameters
        ----------
        inpcrd : str
            File containing coordinates for next round
        states : np array
            Scoring of each macrostate
        n_msm_runs : int array
            List of labels from which generate new inputs
        n_msm_runs_aux : int array
            Number of new runs for each label

        """
        # Determine distribution of new runs according to 'states'
        n_runs = self.n_runs
        print(np.sum(states), len(self.MSM.keep_states),len(states) , states)
        n_msm_runs = np.random.choice(range(len(self.MSM.keep_states)), n_runs, p=states)
        print ('Runs for new epoch:', n_msm_runs)
        
        ## OPTION 1: Use weights to randomly choose a frame and create a corresponding .gro file
        #n_msm_runs_aux = np.zeros(len(self.MSM.keep_states))
        ##print(list(n_msm_runs))#print(self.MSM.keep_states)
        #for i in self.MSM.keep_states:
        #    n_msm_runs_aux[i] = len(analyzer_lib.list_duplicates_of(list(n_msm_runs), i))
        #    #n_msm_runs_aux.append(len(n_msm_runs[ n_msm_runs == i ]))#.count(i)
        #inputs = analyzer_lib.gen_input_weights(n_msm_runs_aux, self.labels_all, self.trajs, tprs)

        # OPTION 2: Use a dict containing all trajs, labels and frames to pick randomly a frame
        state_kv = {}
        inputs = []
        for n in n_msm_runs:
            if n not in state_kv.keys():
                print ("Building entry for microstate %g"%n)
                #analyzer_lib.gen_dict_state(n, self.trajs, state_kv)
                self.gen_dict_state(n, state_kv)
            try:
                traj, frame, which_tr = random.choice(state_kv[n])
                tpr = tprs[which_tr]
                analyzer_lib.map_inputs(traj, n, frame, inputs, tpr)
            except IndexError:
                print('Error in gen_input, check n_msm_runs and MSM.keys')

        return inputs

    def gen_dict_state(self, s, state_kv):
        """
        Generates dictionary entry for state s
    
        Parameters
        ----------
        s : str, int
            The key for the state
        state_kv : dict
            Dictionary containing all trajectories and corresponding frames

        """
        state_kv[s] = []
        n = 0
        for t in self.trajs:
            n +=1
            try:
                ivals = analyzer_lib.list_duplicates_of(t.distraj, self.MSM.keys[s])
                for i in ivals:
                    state_kv[s].append([t.mdt, i, n-1])
            except KeyError:
                pass