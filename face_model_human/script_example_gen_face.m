% run from within face_model_human directory

addpath(['.' filesep 'Model'])

model_data = ['.' filesep 'Model_Data' filesep 'Human_Face_Model_Data.mat'];output_res = [360 250];

model = AAM_Model(model_data); % instance of class AAM_model

%% gen random faces
n_faces = 1;
params = model.gen_random_params(n_faces); % function reads this as gen_random_params(self/model, n_faces)
[im, landmarks] = model.gen_image_param(params, output_res);

%%
for i = 1:n_faces
    figure;imshow(uint8(im(:,:,:,i)))
    
    % plot corresponding landmarks
    if exist('landmarks')
        figure; scatter(landmarks(:, 1, i), -landmarks(:, 2, i));
    end
end


%% messing around with coordinates
setDiskPaths


imPath = [famPath filesep 'Backgroundremoval_SegAny' filesep 'FamFaces_Pt_P87CS_Processed'];
% ims = Utilities.readInFiles(imPath, 'jpg');
ims = Utilities.readInFiles(imPath);

xlsFile = 'imageJ_landmarkcoords_famfaces'; % note these are for aligned, scaled, and BG free images

% keep sheet name same as stim name and put this in loop
n_ims = length(ims);
p_fam = [];

for im_idx = 1:n_ims

% grab the image - needs to be in a cell array for model functions to read
im{1} = imread([ims(im_idx).folder filesep ims(im_idx).name]); % some img I marked landmarks for in imageJ

dashpos = strfind(ims(im_idx).name, '_');
sheet = ims(im_idx).name(1:dashpos-1);

manual_landmarks = xlsread([famPath filesep xlsFile], sheet);


% confirmation 
% figure; scatter(manual_landmarks(:, 1), -manual_landmarks(:, 2));

% produce texture with landmarks marks
m_param = model.compute_param(im, manual_landmarks);
p_fam(im_idx, :) = m_param(1:100);

% % reproduce image with SA params
% [m_im, m_landmarks] =  model.gen_image_param(m_param, output_res);
% figure; imshow(uint8(m_im));
% figure; scatter(m_landmarks(:, 1), -m_landmarks(:, 2));

end

save([famPath filesep 'FeatureMatching' filesep 'params_fam_p87CS_100d.mat'], 'p_fam');


%% renaming fam images

setDiskPaths


imPath = [famPath filesep 'Backgroundremoval_SegAny' filesep 'FamFaces_Pt_P87CS_2'];


files = Utilities.readInFiles(imPath);
startnum = 3000;

changeNum = false;
convertToJPG = true;
for f = 1:numel(files)
    
    fn = [files(f).folder filesep files(f).name];
    fname = files(f).name;
    
    dashpos = strfind(files(f).name, '.');
    suffix = fname(dashpos+1:end);
    prefix = fname(1:dashpos-1);
    % converting to jpgs
    if ~strcmpi(suffix, 'jpg') && convertToJPG
        fname(dashpos+1:end) = 'jpg';
        
        movefile([files(f).folder filesep files(f).name], [files(f).folder filesep fname]);
        
%         im = imread(fn);
%         imwrite(im, [files(f).folder filesep fname]);
        
    end
    
    % change numbers
    if changeNum
        newfname = [num2str(startnum+f) '.' suffix];
        movefile([files(f).folder filesep fname], [files(f).folder filesep newfname]);
    end
end


















