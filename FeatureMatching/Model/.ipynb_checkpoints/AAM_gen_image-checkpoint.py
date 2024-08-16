# -*- coding: utf-8 -*-
import numpy as np
from matplotlib import path
from scipy.interpolate import interpn
from scipy import ndimage

#%%
def get_randn_param(n, data):
    # p = get_randn_param( n, data)
    # generate gaussian distributed random parameters
    # n: number of parameter vectors to generate
    # cut off to the max/min values 
    
    rng = np.random.default_rng()
    ndim = len(data.score_mean)
    p = np.zeros((ndim, n))
    for i in range(ndim):
        mean_1dim = data.score_mean[i]
        std_1dim = data.score_std[i]
        max_1dim = data.score_max[i]
        min_1dim = data.score_min[i]
        
        for j in range(n):
            r = rng.standard_normal() * std_1dim + mean_1dim
            while r > max_1dim or r < min_1dim:
                r = rng.standard_normal(1) * std_1dim + mean_1dim
            
            p[i,j] = r
            
    return p



#%%
def inpolygon(xq, yq, xv, yv):
    # True for points inside a polygonal region
    # xq, yq: points
    # xv, yv: polygon
    
    shape = xq.shape
    n = xq.size
    xq = xq.flatten()
    yq = yq.flatten()
    in_rect = (xq>xv.min()) & (xq<xv.max()) & (yq>yv.min()) & (yq<yv.max())
    q = np.hstack((xq[in_rect,np.newaxis],yq[in_rect,np.newaxis]))
    p = path.Path([(xv[i], yv[i]) for i in range(xv.shape[0])])
    tf = np.zeros(n,dtype=bool)
    tf[in_rect] = p.contains_points(q)
    
    return tf.reshape(shape)

#%%
def polyarea(x,y):
    return 0.5*np.abs(np.dot(x,np.roll(y,1))-np.dot(y,np.roll(x,1)))

#%%

def im_mix(im1, im2, alpha, sigma = 0):
# im_mix(im1, im2, alpha, sigma)
# mix image1 and image2 according to alpha channel
# INPUT: im1, im2, images of same resolution
#             alpha, matrix specify weights of the mix for each pixel
#             sigma: smooth alpha to avoid sharp edge, sigma of Gaussian smooth in pixel

    im1 = im1.astype(np.double)
    im2 = im2.astype(np.double)
    alpha = alpha.astype(np.double)
  
    if sigma >0:
        alpha = ndimage.gaussian_filter(alpha, sigma)
    
    im = im1*alpha + im2*(1-alpha) 
    
    return im

#%%
def im_mask_smooth_edge( im, mask, sigma):
#  make smooth transition along mask edge
#  to gray background 
#  INPUT: 
#      im: image
#      mask: binary array as same size of image 
#                 or polygon defined by sequence of points (n, 2)
#      sigma: sigma of Gaussian smooth in pixel
# 

    y, x = im.shape[0:2]
    
    if sigma is None:
        sigma = y*0.005
    
    bg = np.zeros((1,1,1))+128
    
    if mask.dtype != bool:
        # if it's not bool, it should be a polygon.
        [xx, yy] = np.meshgrid(range(x),range(y))
        mask = inpolygon(xx.flatten(), yy.flaten(), mask[:,0], mask[:,1])
        mask = mask.reshape((y, x))

    
    mask = ndimage.gaussian_filter(mask.astype(np.double), sigma)
    
    mask = (mask - 0.5)*2
    mask[mask<0] = 0
    mask = mask[:,:,np.newaxis]
    
    im = im_mix(im, bg, mask)
    
    return im

#%%
def xy2grid(x, y):
    
    if np.isscalar(x) or len(x)==1 :
        x_grid, y_grid = np.meshgrid(range(x), range(y))
        
    elif len(x.squeeze().shape)==1 :
        # if source coordinate is one dimension
        x_grid, y_grid = np.meshgrid(x, y)
        
    else:
        x_grid = x
        y_grid = y
    
    return x_grid, y_grid

#%%
def affine_transform_triangles( vs, vt, facets ):
    #AFFINE_TRANSFORM_TRIANGLES 
    # vs: vertices source
    # vt: vertices target
    # facets: index of facets
    
    nv = vs.shape[0] 
    vs = np.hstack( (vs, np.ones((nv,1))) )
    vt = np.hstack( (vt, np.ones((nv,1))) )
    
    n = facets.shape[0]
    aft = np.zeros((3,3,n))
    for i in range(n):
        fi = facets[i,:]-1
        tmp = np.linalg.lstsq(vs[fi,:], vt[fi,:], rcond = None)
        aft[:,:,i] = tmp[0]
    return aft

#%%
def draw_texture( source_texture, source_mark, source_x, source_y, \
                 target_mark, target_x, target_y, facets ):
# DRAW_TEXTURE

    # input check
    if len(source_mark.shape)==1 or source_mark.shape[1]==1:
        source_mark = source_mark.reshape((-1,2), order = 'F')
        
    if len(target_mark.shape)==1 or target_mark.shape[1]==1:
        target_mark = target_mark.reshape((-1,2), order = 'F')

    if isinstance(source_texture[0], np.integer):
        source_texture = np.single(source_texture)
    
    #
    sx, sy = xy2grid(source_x, source_y)
    tx, ty = xy2grid(target_x, target_y)
    target_resy, target_resx = tx.shape
    
    # compute affine transformation
    aft = affine_transform_triangles(target_mark, source_mark, facets)
    
    xt = np.full(tx.shape, np.nan) # target x transformed in source coordinate
    yt = np.full(ty.shape, np.nan)

    ##
    nf = facets.shape[0]
    for i in range(nf):
        xv = target_mark[facets[i,:]-1,0]
        yv = target_mark[facets[i,:]-1,1]
        pt_ind = inpolygon(tx, ty, xv, yv) # locate points belong to each facets in target image
        tm = np.hstack( (tx[pt_ind, np.newaxis], ty[pt_ind, np.newaxis], np.ones((pt_ind.sum(),1))) )
        xyt = np.dot(tm, aft[:,:,i]) # piece-wise transform
        xt[pt_ind] = xyt[:,0]
        yt[pt_ind] = xyt[:,1]
    
    ##
    nch = source_texture.shape[2]
    im = np.full((target_resy, target_resx, nch),128)
    
    # fill image
    pt_ind_all = ~np.isnan(xt) & (xt>sx[0,0]) & (xt<sx[0,-1]) & (yt>sy[0,0]) & (yt<sy[-1,0])
    for i in range(nch):
        pti = np.hstack((yt[pt_ind_all,np.newaxis], xt[pt_ind_all,np.newaxis]))
        im[pt_ind_all,i] = interpn((sy[:,0], sx[0,:]), source_texture[:,:,i], pti)
    
    return im, pt_ind_all

#%%
def AAM_gen_image( p_id_mark, p_id_texture, model, output_res ):
    # synthesize_image
    
    # buffer size
    resx = output_res[1]*model.image.upsampling
    resy = output_res[0]*model.image.upsampling
    
    ## mark
    nmark = len(p_id_mark)
    neutral_mark = model.id_mark.mean.T + np.dot(model.id_mark.eig_vector[:,:nmark], p_id_mark)
    neutral_mark = neutral_mark.reshape((-1,2), order='F')
    
    # scale face area to proportion of images
    peri_ind = np.hstack((model.mark_group.rim, model.mark_group.rim[0]))-1
    area = polyarea(neutral_mark[peri_ind,0], neutral_mark[peri_ind,1])
    scaling_factor = np.sqrt(resx * resy * 0.6 / area)
    
    dy = model.image.y_offset_factor * resx
    
    neutral_mark[:,0] = neutral_mark[:,0] * scaling_factor + resx/2
    neutral_mark[:,1] = neutral_mark[:,1] * scaling_factor + resy/2 + dy
    
    ## texture
    ntexture = len(p_id_texture)
    resy_texture, resx_texture = model.x_texture.shape
    neutral_texture = model.id_texture.mean.T + np.dot(model.id_texture.eig_vector[:,:ntexture], p_id_texture)
    neutral_texture = neutral_texture.reshape((resy_texture, resx_texture, 3), order='F')
    
    im, mask = draw_texture( neutral_texture, model.id_mark.mean, model.x_texture, model.y_texture, \
                      neutral_mark, resx, resy, model.facets )
    
    ## smooth edge
    if model.image.smooth_edge:
        sigma = model.image.smooth_factor*resx
        # edge = neutral_mark[model.mark_group.rim,:]
        im = im_mask_smooth_edge(im, mask, sigma)
    
    ## more efficient, lower quality of upsampling
    # if model.image.upsampling>1:
    #     im = ndimage.zoom(im, model.image.upsampling)
    
    
    landmark = neutral_mark
    
    return im, landmark