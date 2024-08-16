# -*- coding: utf-8 -*-

import pickle
import os
import numpy as np

from scipy import stats


#%%
p_fam = np.load('./all_face_features_60_20d_norm.npy')

#%% get best faces

import fnmatch

# data_folder = 'results_eps001_dev1'; deviation = 1
data_folder = 'results_eps001_dev15'; deviation = 1.5

n_repeats = 5


def get_costs(filename):
    
    i = int(filename[6:8])
    j = int(filename[9:11])
    
    file = os.path.join(data_folder,filename)
    
    with open(file, 'rb') as f:
        r = pickle.load(f)

    cost = r.fun
    
    return i, j, cost


files = os.listdir(data_folder)
files = fnmatch.filter(files,'*.pickle')

nfile = len(files)

costs = np.zeros((60, n_repeats)) 
for k in range(nfile):
    i, j, c = get_costs(files[k])
    costs[i,j] = c

index_selected = np.argmin(costs,axis=1)

#%% get results

d = np.zeros((60,20))

for i in range(60):
    
    j = index_selected[i]
    filename = f'result{i:02d}_{j:02d}.pickle'
    file = os.path.join(data_folder,filename)
    
    with open(file, 'rb') as f:
        r = pickle.load(f)
    
    x = r.x
    print(np.std(x))
    x = x/np.sqrt(np.inner(x, x)/len(x))*deviation
    #x = x.reshape((1,-1))
    d[i,:]=x
    
np.save('all_similar_face_pair_deviation_60_20d.npy',d)