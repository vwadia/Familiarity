# -*- coding: utf-8 -*-
# gen face similar to target face (fixed small distance in face space) 
# and low level feature matched

# Distribution difference measure: KL divergence
# Match outer contour also

import time
t_start = time.time()


#%% 

import numpy as np
# import concurrent.futures
# import matplotlib.pyplot as plt

from scipy import optimize

from skimage.color import rgb2hsv, rgb2gray
import pickle


#%% input
import sys
if len(sys.argv)>1:
    iface = int(sys.argv[1])
else:
    iface = 0

#%% check where the code is running
import os
if os.name=='nt':
    sys.path.append('D:/github/AAM')
    model_data = 'D:/Project2 memory/2D_monkey_face_model/Model_Data/Face_Model_Data.mat';
else:
    model_data = '../data/Face_Model_Data.mat'


#%% optimization options

# op_method = 'BFGS'
# options_optimize = {'disp': True, 'maxiter': 1000, 'gtol': 1e-03, 'eps': 1e-02}

op_method = 'Nelder-Mead'
options_optimize = {'disp': True, 'maxiter': 1000}


#%% initialize face generator
from Model.AAM_Model import AAM_Model

model = AAM_Model(model_data)
output_res = [300, 300]
options_face_gen = {'normalized': True, 'ndim_shape': 10}

#%% familiar faces params (shape-appearance feature)
p_fam = np.load('all_face_features_60_20d_norm.npy')
n_fam = p_fam.shape[0]
n_dim = p_fam.shape[1]

#%%
rng = np.random.default_rng()


#%% get image features
bin_hue = np.arange(0,1.01,0.01)
bin_low = np.arange(-10,10,0.5)

def get_image_features_distribution(im):
    
    im_gray = rgb2gray(im.squeeze())
    
    lum = im_gray.mean()
    contrast = np.square(im_gray-128).mean()
    
    im_hsv = rgb2hsv(im.squeeze())
    hue = im_hsv[:,:,0].flatten()
    hue_pdf = np.histogram(hue, bin_hue)[0]/len(hue)

    return lum, contrast, hue_pdf


#%% optimize



def find_low_level_matched_face(p_fam1, deviation):
    #familiar face low level features
    
    im_fam, landmarks_fam = model.gen_image_param(p_fam1.reshape((1,-1)), output_res, options_face_gen)
    rim_fam = landmarks_fam[model.data.mark_group.rim,:]
    
    lum_fam, contrast_fam, hue_fam_pdf = get_image_features_distribution(im_fam)
    hue_fam_entropy = - np.dot(hue_fam_pdf, np.log2(hue_fam_pdf+1e-6))
    
    def cost_fun(dp):
        
        ## constraint in distance of deviation vector
        dp = dp/np.sqrt(np.inner(dp, dp)/len(dp))*deviation
        
        p = p_fam1 + dp
                
        im, landmarks = model.gen_image_param(p.reshape((1,-1)), output_res, options_face_gen)
        rim = landmarks[model.data.mark_group.rim,:]
        
        lum, contrast, hue_pdf = get_image_features_distribution(im)
        
        # contrast
        cost1 = abs(contrast_fam - contrast)/100
        
        # luminance
        cost2 = abs(lum_fam - lum)/10
        
        # distribution hue
        cost3 = (-np.dot(hue_fam_pdf, np.log2(hue_pdf+1e-6)) - hue_fam_entropy)*10
        
        # rim
        cost4 = np.mean((rim_fam - rim)**2)/200
        
        # total cost
        cost = cost1 + cost2 + cost3 + cost4
        
        # print(f'cost:{cost:.4f},contr:{cost1:.4f},lum:{cost2:.4f},hue:{cost3:.4f},rim:{cost4:.4f}')
        # print(time.time())
        
        return cost
    
    dp0 = rng.normal(scale = deviation, size = p_fam1.shape)    
    
    print ('start optimization...')
    res = optimize.minimize(cost_fun, dp0, method = op_method, options = options_optimize)
    
    return res


#%%
deviation = 1.5

# results = []
for i in range(5):
    
    print (f'face{iface}:')
    p_fam1 = p_fam[iface,:]
    r = find_low_level_matched_face(p_fam1, deviation)
    # results.append(r)
    
    with open(f'./results/result{iface:02d}_{i:02d}.pickle', 'wb') as f:
        pickle.dump(r, f)


t_end = time.time()
t_diff = t_end - t_start
print(f'elapsed time {t_diff} seconds')