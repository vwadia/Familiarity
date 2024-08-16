% script to readin marked points and produce reconstructed faces from the
% models


% set paths
run G:\SUAnalysis\setDiskPaths

% initialize SA model
modelPath = [famPath filesep 'face_model_human' filesep 'Model'];
addpath(modelPath)

model_data = [famPath filesep 'face_model_human' filesep 'Model_Data' filesep 'Human_Face_Model_Data.mat'];
output_res = [360 250];

model = AAM_Model(model_data); % instance of class AAM_model


%% read in coords

useAligned = false;

if useAligned
    
%     imPath = [famPath filesep 'BackgroundRemoval_SegAny' filesep 'SingleIm_Aligned'];
%     imPath = [famPath filesep 'BackgroundRemoval_SegAny' filesep 'CelebFaces_pt_Varun_Aligned_10'];
    imPath = [famPath filesep 'BackgroundRemoval_SegAny' filesep 'P86CS' filesep 'FamFaces_Pt_P86CS_Aligned'];
    ims = Utilities.readInFiles(imPath);
    
%     ptsPath = [famPath filesep 'BackgroundRemoval_SegAny' filesep 'SingleIm_markedPts_Aligned'];
%     ptsPath = [famPath filesep 'BackgroundRemoval_SegAny' filesep 'markedPts_pt_Varun_Aligned_10'];
    ptsPath = [famPath filesep 'BackgroundRemoval_SegAny' filesep 'P86CS' filesep 'markedPts_pt_P86CS'];
    m_pts = Utilities.readInFiles(ptsPath, 'mat');   
    
    % for SA images
    outPath = [famPath filesep 'BackgroundRemoval_SegAny' filesep 'P86CS' filesep 'SA_Aligned'];
else
    
%     imPath = [famPath filesep 'BackgroundRemoval_SegAny' filesep 'SingleIm'];
%     imPath = [famPath filesep 'BackgroundRemoval_SegAny' filesep 'CelebFaces_pt_Varun_Processed_10'];
%     imPath = [famPath filesep 'BackgroundRemoval_SegAny' filesep 'FamFaces_pt_P87CS_Processed'];
%     imPath = [famPath filesep 'BackgroundRemoval_SegAny' filesep 'FamFaces_pt_P87CS_2'];
    % imPath = [famPath filesep 'BackgroundRemoval_SegAny' filesep 'P92CS' filesep 'FamFaces_Pt_P92CS'];
    imPath = [famPath filesep 'BackgroundRemoval_SegAny' filesep 'P98CS' filesep 'FamFaces_P98CS_Raw'];
    ims = Utilities.readInFiles(imPath, 'jpg');
    
%     ptsPath = [famPath filesep 'BackgroundRemoval_SegAny' filesep 'SingleIm_markedPts'];
%     ptsPath = [famPath filesep 'BackgroundRemoval_SegAny' filesep 'markedPts_pt_Varun_Processed_10'];
%     ptsPath = [famPath filesep 'BackgroundRemoval_SegAny' filesep 'markedPts_pt_P87CS_Processed'];
%     ptsPath = [famPath filesep 'BackgroundRemoval_SegAny' filesep 'markedPts_pt_P87CS_2'];
    % ptsPath = [famPath filesep 'BackgroundRemoval_SegAny' filesep 'P92CS' filesep 'markedPts_Pt_P92CS'];
    ptsPath = [famPath filesep 'BackgroundRemoval_SegAny' filesep 'P98CS' filesep 'markedPts_Pt_P98CS'];
    m_pts = Utilities.readInFiles(ptsPath, 'mat'); 
%     m_pts = m_pts(4);
    
    % for SA images
%     outPath = [famPath filesep 'BackgroundRemoval_SegAny' filesep 'SA_P87CS_2'];
    % outPath = [famPath filesep 'BackgroundRemoval_SegAny' filesep 'SA_P92CS'];
    outPath = [famPath filesep 'BackgroundRemoval_SegAny' filesep 'SA_P98CS'];
end

if ~exist(outPath, 'dir')
    mkdir(outPath)
end

%%
% perform interpolation to get final 12 points
toSave = true;
display = false;
for i = 1:length(m_pts)

    [~, n1, ext1] = fileparts([m_pts(i).folder filesep m_pts(i).name]);
    [~, n2, ext2] = fileparts([ims(i).folder filesep ims(i).name]);
    assert(strcmp(n1, n2));
    
    % interpolate
    file80 = m_pts(i);
    pts = get_landmark92(file80, toSave); % point saving happens here
    
    [loaded_im, cmap] = imread([ims(i).folder filesep ims(i).name]); % some img I marked landmarks for in imageJ
    
    if ~isempty(cmap)
        loaded_im = ind2rgb(loaded_im, cmap);
    end
%     if size(loaded_im, 1) ~= output_res(1) && size(loaded_im, 2) ~= output_res(2)
%         im{1} = imresize(loaded_im, output_res);
%     else
        im{1} = loaded_im;
%     end
    % produce new face
    m_param = model.compute_param(im, pts);

    % reproduce image with SA params
    [m_im, m_landmarks] =  model.gen_image_param(m_param, output_res);
    m_im = uint8(m_im);
    if display
        f = figure; imshow(m_im);
    end
%     keyboard
    p_fam(i, :) = m_param(1:100);
    imwrite(m_im, [outPath filesep [n2 ext2]])
end

if toSave
    save([famPath filesep 'FeatureMatching' filesep 'params_fam_P98_100d.mat'], 'p_fam');
    % save([famPath filesep 'FeatureMatching' filesep 'params_fam_P92_100d.mat'], 'p_fam');
end

