% familiar face play 
% seeing which algorithm produced the best feature matched face for the first 
% familiar face (aaron paul)

% vwadia March 2023

%% paths
setDiskPaths

res_path = [famPath filesep 'FeatureMatching'];
addpath([famPath filesep 'face_model_human' filesep 'Model'])

model_data = [famPath filesep 'face_model_human' filesep 'Model_Data' filesep 'Human_Face_Model_Data.mat'];

model = AAM_Model(model_data); % instance of class AAM_model



nm_results  = load([res_path filesep 'best_faceparams_AaronPaul_NM.mat']); 
% bfgs_results  = load([res_path filesep 'best_faceparams_AaronPaul_BFGS.mat']); 
bfgs_results  = load([res_path filesep 'best_faceparams_AaronPaul_1k_BFGS.mat']); 

load([res_path filesep 'params_fam_100d.mat'])
target = p_fam(1, :);

%%
     
% NM
% nm_dist = vecnorm(nm_results.d' - repmat(target, [size(nm_results.d, 1) 1])');
% nm_dist = squareform(pdist([target; nm_results.d], 'cosine'));
% nm_dist = nm_dist(2:end, 1);
% [nm_dist_min, nm_idx] = min(nm_dist);

% BFGS
% bfgs_dist = vecnorm(bfgs_results.d' - repmat(target, [size(bfgs_results.d, 1) 1])');
bfgs_dist = squareform(pdist([target; bfgs_results.d], 'cosine'));
bfgs_dist = bfgs_dist(2:end, 1);
[bfgs_dist_min, bfgs_idx] = min(bfgs_dist);


%% produce faces with those values

output_res = [360 256];
options.normalized = true;

for i = 1:size(bfgs_results.d, 1)
    
%     params = nm_results.d(i, :);
    params = bfgs_results.d(i, :);
    [im, landmarks] = model.gen_image_param(params, output_res, options);
    figure; imshow(uint8(im)); title(['Face ' num2str(i)]);

end

%% actual face

output_res = [360 256];
for i = 1%:14
    target = p_fam(1, :);
    options.normalized = false;
    params = target;
    [im, landmarks] = model.gen_image_param(params, output_res, options);
    figure; imshow(uint8(im)); title(['Fam Face ' num2str(i)]);
%     
end

%%
imPath = [famPath filesep 'FamFacesProcessed']; 
ims = Utilities.readInFiles(imPath, 'jpg');

for i = 15:length(ims)
    
%     strpos = strfind(ims(i).name, '_');
%     name = ims(i).name(1:strpos-1);
%         newname = sprintf('%03d', str2num(i+99));

    newname = sprintf('%03d', i+100);
    
    
    movefile([ims(i).folder filesep ims(i).name], [ims(i).folder filesep [newname '.jpg']]);
    
    
end

