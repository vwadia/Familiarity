# -*- coding: utf-8 -*-

#%% 
import numpy as np
import concurrent.futures
# import matplotlib.pyplot as plt

import torch
from torchvision import transforms, models

import scipy.io as sio
from scipy import stats
from scipy import optimize
from scipy.spatial.distance import pdist

#%% optimization options

options_optimize = {'disp': True, 'maxiter': 100, 'gtol': 1e-03, 'eps': 1e-03}
workers = 32

#%% initialize face generator
from Model.AAM_Model import AAM_Model
model_data = '../data/Face_Model_Data.mat'
model = AAM_Model(model_data)
output_res = [300, 300]
options_face_gen = {'normalized': True}

#%% initialize NN
base_model_name = "alexnet"
feature_extract = True

model_alexnet = models.alexnet(pretrained=True)
for param in model_alexnet.parameters():
    param.requires_grad = False
model_alexnet.eval()

#%% familiar faces params (shape-appearance feature)
p_fam = np.load('familiar_face_features_36_norm.npy')
n_fam = p_fam.shape[0]
n_dim = p_fam.shape[1]

output_res_n = [output_res]*n_fam
options_face_gen_n = [options_face_gen]*n_fam

d_fam = pdist(p_fam)
var_fam = np.var(p_fam)

#%% 
def im_preproces(im):
    
    im = im.astype(np.single)
    im[im>255] = 255
    im[im<0] = 0
    
    im = np.moveaxis(im, (-1,-2), (0,1))
    
    mean = np.array([0.485, 0.456, 0.406], dtype=np.single).reshape((1,3,1,1))
    std = np.array([0.229, 0.224, 0.225], dtype=np.single).reshape((1,3,1,1))
    im = (im/255.0-mean)/std
    im = torch.from_numpy(im)
    return im
    
#%% familiar face low level features

im_fam = model.gen_image_param(p_fam, output_res, options_face_gen)
im_fam = im_preproces(im_fam)
features_low_fam = model_alexnet.features[0](im_fam)
features_low_fam = features_low_fam.data.numpy().reshape((n_fam,-1))
d_fam_low = pdist(features_low_fam)


#%%
def cost_fun(p):
    
    p = p.reshape((n_fam,n_dim))
    
    params_n = [p[np.newaxis,i,:] for i in range(n_fam)]
    with concurrent.futures.ThreadPoolExecutor(max_workers=workers) as executor:
        im = np.concatenate(list(executor.map(model.gen_image_param, params_n, output_res_n, options_face_gen_n)), axis=3)

    im = im_preproces(im)
    features_low = model_alexnet.features[0](im)
    features_low = features_low.reshape((n_fam,-1))
    
    ## shape-appearance feature
    # variance
    var = np.var(p);
    cost1 = np.mean(np.square(var-var_fam)) # mse
    cost2 = abs(np.mean(var) - np.mean(var_fam)) # total variance
    
    # distribution each dimension
    ps = np.zeros((n_dim + 3));
    for i in range(n_dim):
        ksresult = stats.kstest(p[:,i],p_fam[:,i]);
        ps[i] = ksresult.pvalue
    
    # pairwise distance distribution
    d = pdist(p);
    ksresult = stats.kstest(d,d_fam);
    ps[n_dim] = ksresult.pvalue
    
    ## low level features
    # pairwise distance distribution
    d_low = pdist(features_low);
    ksresult = stats.kstest(d_low, d_fam_low);
    ps[n_dim+1] = ksresult.pvalue
    
    # distribution low level feature
    ksresult = stats.kstest(features_low.flatten(), features_low_fam.flatten());
    ps[n_dim+2] = ksresult.pvalue
    
    # compute cost3 (p value for all distributions)
    p_min = ps.min()
    
    if p_min < 0.05:
        cost3 = 1/(p_min+0.00001)
    else:
        cost3 = 0
    
    # total cost
    cost = cost1 + cost2 + cost3
    
    print(f'cost:{cost:.4f},mse:{cost1:.4f},var:{cost2:.4f},pd_low:{ps[n_dim+1]:.4f},p_low:{ps[n_dim+2]:.4f}')
    # print(f'cost:{cost:.4f}')
          
    return cost

#%% optimize
rng = np.random.default_rng()

p0 = np.zeros(p_fam.shape)
for i in range(n_dim):
    p0[:,i] = p_fam[rng.permutation(n_fam),i]

res = optimize.minimize(cost_fun, p0, method = 'BFGS', options = options_optimize)

