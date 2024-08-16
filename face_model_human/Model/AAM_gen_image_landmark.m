function  im = AAM_gen_image_landmark( landmark, p_id_texture, model, output_res )
% AAM gen image by input of actual facial landmarks
%

% buffer size
resx = output_res(2)*model.image.upsampling; 
resy = output_res(1)*model.image.upsampling;


%% mark
neutral_mark = landmark;

%% texture
[resy_texture, resx_texture] = size(model.x_texture);
neutral_texture = model.id_texture.mean' + model.id_texture.eig_vector * p_id_texture ;
neutral_texture = reshape(neutral_texture, resy_texture, resx_texture, 3);

im = draw_texture( neutral_texture, model.id_mark.mean, model.x_texture, model.y_texture, neutral_mark, resx, resy, model.facets );

%% smooth edge
if model.image.smooth_edge
    sigma = model.image.smooth_factor*resx;
    edge = neutral_mark(model.mark_group.rim,:);
    im = im_mask_smooth_edge(im, edge, sigma);
end

%%
if model.image.upsampling>1
    im = imresize(im, output_res);
end

end

