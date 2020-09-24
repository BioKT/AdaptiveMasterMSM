"""
This file is part of the AdaptiveMasterMSM package.

"""

import numpy as np
import random
#import mdtraj as md

def gen_input_weights(n_msm_runs, labels_all, trajs):
    """
        
    Parameters
    ----------
    n_msm_runs : int array
        Number of runs for each macrostate
    labels_all : list
        List of all discretized trajectories
    trajs : List
        List of instances of TimeSeries for all input trajectories

    """

    counter = np.zeros(len(n_msm_runs))
    # generar weights
    w = np.zeros((len(labels_all), len(n_msm_runs)))
    for j, x_j in enumerate(labels_all):
        for i in range(len(n_msm_runs)):
            # row: which traj, column: which label
            index = list_duplicates_of(list(x_j), i)
            w[j][i] = len(index)
        if np.sum(w[j]) > 0.0: w[j] /= np.sum(w[j])

    for i in range(len(n_msm_runs)):
        while True:
            which_tr_index = random.choices(range(len(labels_all)), weights=w[:,i])
            which_tr = labels_all[which_tr_index[0]]
            # obtain positions in traj where label 'i' is found
            index = list_duplicates_of(list(which_tr), i)
            if len(index) is 0: continue # esto es redundante si weights funciona
            j = random.choice(index)
            map_inputs(trajs, which_tr_index[0], j)
            counter[which_tr[j]] += 1
            if n_msm_runs[i] <= counter[which_tr[j]]: break
        
    print(np.sum(counter), counter)

def map_inputs(trajs, which_tr_index, state):
    """
        
    Parameters
    ----------
    which_tr_index : int
        Which trajectory identified by index
    state : int
        Index corresponding to the initial state for the new input

    """

    #print(type(self.trajs[which_tr_index].mdt), type(self.trajs[which_tr_index].mdt[state]))
    traj = trajs[which_tr_index].mdt
    fn = '%s_%s.gro'%(which_tr_index,state)
    traj[state].save_gro(fn)


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
