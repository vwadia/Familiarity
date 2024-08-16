import numpy as np
import scipy.io as sio


from .AAM_gen_image import AAM_gen_image, get_randn_param

#%% AAM_Model
class AAM_Model():
    # AAM_MODEL Active Apearance Model
    
    def __init__(self, data_file):
        model_data = sio.loadmat(data_file, struct_as_record=False, squeeze_me=True)
        self.data = model_data['model'];
        self.npc_mark = model_data['npc_id_mark'] # number of PC
        self.npc_texture = model_data['npc_id_texture']
        
    #%% helper
    
    def gen_random_params(self, n):
        # generate random paramters
        # n: number of samples to generate
        p_id_mark = get_randn_param( n, self.data.id_mark)
        p_id_texture = get_randn_param( n, self.data.id_texture)
        return np.hstack((p_id_mark.T, p_id_texture.T))

        
    #%% gen_image_param ()
    # generate image from provided parameters
    def gen_image_param(self, params, output_res, options=None):
    # [im_syn, landmarks] = gen_image_param(self, params, output_res, options = None)
    # INPUT:
    #    params: shape-appearance paramters 
    #    output_res: output resolution
    # 
    #    [optional] 
    #        options['ndim_shape']: number of shape dimension
    #        options['normalized']: true or false, params is normalized or not
        
        if isinstance(options, dict):
            nmark = options.get('ndim_shape', self.npc_mark)
            normalized_params = options.get('normalized', False)
        else:
            nmark = self.npc_mark
            normalized_params = False
        
        ntexture = params.shape[1] - nmark
        
        # if params are normalized, recover its scale
        if normalized_params:
            p_std = np.hstack((self.data.id_mark.score_std[0:nmark], self.data.id_texture.score_std[0:ntexture])).reshape((1,-1))
            p_mean = np.hstack((self.data.id_mark.score_mean[0:nmark], self.data.id_texture.score_mean[0:ntexture])).reshape((1,-1))
            params = params * p_std + p_mean
        
        n = params.shape[0]
        im_syn = np.zeros(output_res + [3, n]);
        landmarks = np.zeros((len(self.data.id_mark.mean)//2, 2, n));
                
        for i in range(n):
            #print(f'Generating images {i}/{n}')
            p_id_mark = params[i,:nmark]
            p_id_texture = params[i,nmark:]
            im_syn[:,:,:,i], landmarks[:,:,i] = AAM_gen_image( p_id_mark, p_id_texture, self.data, output_res)
        
        im_syn[im_syn>255] = 255
        im_syn[im_syn<0] = 0
        
        return im_syn, landmarks