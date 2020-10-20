"""
This file is part of the AdaptiveMasterMSM package.

"""

import sys, os
import numpy as np
import random
from adaptivemastermsm.launcher import launcher_lib
#import mdtraj as md

def gen_input_weights(n_msm_runs, labels_all, trajs, tprs):
    """
        
    Parameters
    ----------
    n_msm_runs : int array
        Number of runs for each macrostate
    labels_all : list
        List of all discretized trajectories
    trajs : List
        List of instances of TimeSeries for all input trajectories
    tprs : List
        tpr files corresponding to 'trajs' trajectories

    """

    # generar weights
    w = np.zeros((len(labels_all), len(n_msm_runs)))
    for j, x_j in enumerate(labels_all):
        for i in range(len(n_msm_runs)):
            # row: which traj (or list of labels), column: which label
            index = list_duplicates_of(list(x_j), i)
            w[j][i] = len(index)
        if np.sum(w[j]) > 0.0: w[j] /= np.sum(w[j])

    inputs = []
    counter = np.zeros(len(n_msm_runs))
    for i in range(len(n_msm_runs)):
        while True:
            which_tr_index = random.choices(range(len(labels_all)), weights=w[:,i])
            which_tr = labels_all[which_tr_index[0]]
            # obtain positions in traj where label 'i' is found
            index = list_duplicates_of(list(which_tr), i)
            print('ionix',i)
            print(w[:,i])
            if len(index) is 0: break # redundant if 'weigths' works
            j = random.choice(index)
            map_inputs(trajs[which_tr_index[0]].mdt, i, j, inputs, tprs[which_tr_index[0]])
            counter[which_tr[j]] += 1
            if n_msm_runs[i] <= counter[which_tr[j]]: break
        
    #print(np.sum(counter), counter)
    #print(inputs)

    return inputs

def map_inputs(traj, label, frame, inputs, tpr):
    """
        
    Parameters
    ----------
    traj : object
        mdtraj Trajectory object (see in TimeSeries of MasterMSM/traj.py)
    label : int
        Label identifying the initial state of the new input
    frame : int
        Frame corresponding to the initial state for the new input

    """

    fn = '%s_%s.gro'%(label,frame)
    inputs.append(fn)
    #traj[frame].save_gro(fn) #traj[frame].center_coordinates()
    ###tr = 'data/tr_%s.h5'%(frame)
    ###traj.save_hdf5(tr)
    tr = 'data/%s_%s.xtc'%(label,frame)
    traj.save_xtc(tr)
    cmd = "gmx trjconv -s %s -f %s -o %s -dump %s <<EOF\n1\nEOF"%(tpr, tr, fn, traj.time[frame])
    os.system(cmd)

def list_duplicates_of(seq,item):
    start_at = -1
    locs = []
    while True:
        try:
            loc = seq.index(item,start_at+1)
        except ValueError:
            break
        else:
            locs.append(loc)
            start_at = loc
    return locs

