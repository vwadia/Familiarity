% 

clear
% f = function_handle('G:\SUAnalysis\setDiskPaths.m');
% f()
%parpool(32)


[~, host] = system('hostname');
if strcmp(host(1:end-1), 'DWA644201')
    atCedars = 1;
    famPath = 'C:\Users\wadiav\Documents\PYTHON\Familiarity';
    diskPath = 'G:\SUAnalysis';
    
elseif strcmp(host(1:end-1), 'DESKTOP-LJHLIED')
    atCedars = 0;
    famPath = 'C:\Users\varunwadia\Documents\PYTHON\Familiarity';
    diskPath = 'G:\SUAnalysis';
    
else % mac
    atCedars = 0;
    famPath = '/Users/varunwadia/Documents/Familiarity';
    diskPath = '/Volumes/T7/SUAnalysis';
    
end

matchPath = [famPath filesep 'FeatureMatching'];

% patientID = 'P87CS';
% patientID = 'P86CS';
% patientID = 'P92CS';
% patientID = 'P98CS';
% patientID = 'P99CS';
patientID = 'P103CS';

outPath = [matchPath filesep patientID];
if ~exist(outPath, 'dir')
    mkdir(outPath)
end

% arbitrary choice
% dim_ind = [1:10 21:50];
% dim_ind = [1:60]; 
dim_ind = [1:100];


% dim_ind = [1:10 21:30]; % his fam SA features had only 10 shape and 10 app dims

forRecall = true; 

%% faces
load([famPath filesep 'FeatureMatching' filesep 'params_1k_100d.mat'])
% load([famPath filesep 'FeatureMatching' filesep 'params_fam_100d.mat'])
% load([famPath filesep 'FeatureMatching' filesep 'params_fam_P87_100d.mat'])
% load([famPath filesep 'FeatureMatching' filesep 'params_fam_P87_2_100d.mat'])
% load([famPath filesep 'FeatureMatching' filesep 'params_fam_P86_100d.mat'])
load([famPath filesep 'FeatureMatching' filesep 'params_fam_' patientID(1:end-2) '_100d.mat'])

unfam_ind = 1:1000;
allfam_ind = 1001:(1000+size(p_fam, 1));

params = [params(unfam_ind, :); p_fam];


params = bsxfun(@rdivide, params, std(params)); % normalize features

subsampleims = false;
if strcmp(patientID, 'P98CS')
    recallIms = [13 15 18 24 25 38 40 45 47 58 64 83]; % chosen manually per patient
elseif strcmp(patientID, 'P99CS')
    recallIms = [1 2 10 27 38 39 75 87 92 93 98 99]; % these are wrt the fam_stim for that patient not the big set
elseif strcmp(patientID, 'P103CS')
    recallIms = [1 3 12 16 21 23 27 45 50 67 71 93]; % these are wrt the fam_stim for that patient not the big set
end
if exist('recallIms', 'var') && forRecall
    refineSearch = true; % can use to find matches for recall Ims - vwadia July 2024
else
    refineSearch = false;
end

%% sub sample params
if subsampleims || refineSearch
    if ~refineSearch
        sub_sample = sort(randi(length(unfam_ind), [1 500])); % NEED TO FUCKING SAVE THESE INDICES vwadia May2024
    else
        % load('match_by_grad_Nfam60_Nunfam184_1000.mat')
        load([outPath filesep 'match_by_grad_Nfam100_Nunfam120_1000.mat'])
        p_fam = p_fam(recallIms, :);
        allfam_ind = recallIms+length(unfam_ind);
        sub_sample = ind_unfam_sub;
        clearvars ind_unfam_sub ind_fam_sub
%        subsample = 1:250;
    end
    unfam_ind = sub_sample';
    % allfam_ind = unfam_ind+1:(unfam_ind+size(p_fam, 1));
    
    params = params([unfam_ind allfam_ind], :);
    allfam_ind = length(unfam_ind)+1:length(unfam_ind)+size(p_fam, 1);
end
%%

params = params(:, dim_ind);
n_stim = size(params,1);

x_fam = false(n_stim,1); x_fam(allfam_ind) = true;

x_unfam = false(n_stim,1); x_unfam(1:length(unfam_ind)) = true;

%%

% Liang said just make this a large number.
% If left alone the initial cost put out by cost_fun is ~1020
% so made it bigger than that
cost = 2e3; 

for k = size(p_fam, 1)%[33 34 35] %[23 24 25 26] %  number of fam faces you want matched
    
    tic
    n_fam_min = k;

x = true(n_stim,1);
n_fam = length(allfam_ind);
n_unfam = length(unfam_ind);


n = 2000; %number of steps

grad = zeros(n_stim,1);
for i = 1:n
%      tic
 
    % estimate grad
    parfor j=1:n_stim 
%     for j=1:n_stim
        xd = x;
        xd(j) = ~xd(j);
        grad(j) = cost_fun( params, xd, x_fam, x_unfam ) - cost;
    end
    
    if n_fam<=n_fam_min
        grad(x&x_fam) = inf; % finds the indices where both x and x_fam = 1
    end
    [grad_min, ind] = min(grad);
    
    if grad_min<0 && n_unfam>(1.2*n_fam) % give a lil extra for repetitions
        x(ind) = ~x(ind);
        cost = cost + grad(ind);
        
        n_fam = sum(x&x_fam);
        n_unfam = sum(x&x_unfam);
        fprintf('[%d] cost %.5f, Nfam%d, Nunfam%d,\n', i, cost, n_fam, n_unfam); % ~10s per iteration
    
%     elseif n_unfam<=n_fam
%         
%         keyboard
%         break;
        
    else
        break;
    end
%     toc

end
toc

ind_unfam_sub = find(x&x_unfam);
ind_fam_sub = find(x&x_fam);
if exist("recallIms", 'var') && forRecall
    output_file_name = sprintf('ForRecall_match_by_grad_Nfam%d_Nunfam%d_%d.mat', n_fam, n_unfam, length(unfam_ind));
else
    output_file_name = sprintf('match_by_grad_Nfam%d_Nunfam%d_%d.mat', n_fam, n_unfam, length(unfam_ind));
end
save([outPath filesep output_file_name], 'ind_unfam_sub', 'ind_fam_sub');
 
end

%% write out the relevant images 
total_unfam_imgs = 1000;

if exist('sub_sample', 'var') 
    if exist("recallIms", 'var')
        total_unfam_imgs = length(sub_sample);
    else
        total_unfam_imgs = 1000 - length(sub_sample);
    end
    % these are always at the end
    famIm_inds = ind_fam_sub + total_unfam_imgs;
    
    % actually important indices
    unfamIm_inds = unfam_ind(ind_unfam_sub);

else
    
    % these are always at the end
    famIm_inds = ind_fam_sub;
    
    % actually important indices    
    unfamIm_inds = unfam_ind(ind_unfam_sub);
end

% unfamiliar images
imagePathUnfam = [diskPath filesep 'ObjectSpace' filesep 'FamiliarFaceCode' filesep '1000_usedInModel'];
all_unfamImgs = Utilities.readInFiles(imagePathUnfam, 'png');

% familiar images
% imagePathFam = [famPath filesep 'BackgroundRemoval_SegAny' filesep 'FamFaces_Pt_P87CS_2_Processed'];
% imagePathFam = [famPath filesep 'BackgroundRemoval_SegAny' filesep patientID filesep 'FamFaces_Pt_P86CS_Processed'];
% imagePathFam = [famPath filesep 'BackgroundRemoval_SegAny' filesep patientID filesep 'FamFaces_Pt_P92CS_Processed'];
imagePathFam = [famPath filesep 'BackgroundRemoval_SegAny' filesep patientID filesep 'FamFaces_Pt_' patientID '_Processed'];
all_famImgs = Utilities.readInFiles(imagePathFam, 'jpg');

% set up output directory

outPath = [diskPath filesep 'Fam_Task' filesep patientID filesep output_file_name(1:end-4)];

if ~exist(outPath)
    mkdir(outPath)
end
%%
% copy over the images - unfamiliar first
unfamImgs = all_unfamImgs(unfamIm_inds);
startnum = 5000;
for uI = 1:length(unfamImgs)
    fnum = startnum+uI;
    newfnum = str2num(unfamImgs(uI).name(1:end-4))+startnum;
    fname = unfamImgs(uI).name;
    dotpos = strfind(unfamImgs(uI).name, '.');
    suffix = fname(dotpos+1:end);
    if ~exist('recallIms', 'var')
        copyfile([unfamImgs(uI).folder filesep unfamImgs(uI).name], [outPath filesep num2str(newfnum) '.' suffix]);
        % copyfile([unfamImgs(uI).folder filesep unfamImgs(uI).name], [outPath filesep num2str(fnum) '.' suffix]);
    else
        copyfile([unfamImgs(uI).folder filesep unfamImgs(uI).name], [outPath filesep num2str(newfnum) '.' suffix]); % keep the names the same
        % copyfile([unfamImgs(uI).folder filesep unfamImgs(uI).name], [outPath filesep unfamImgs(uI).name]); % keep the names the same

    end
    
end
%%
if ~forRecall
% if ~exist('recallIms', 'var')
    % familiar second
    useoffset = false; % inthis big set these are already offset
    famImgs = all_famImgs(famIm_inds - total_unfam_imgs);

    for fI = 1:length(famImgs)
        fn = famImgs(fI).name;
        if useoffset
            fInd = 3000 + str2num(fn(1:end-4));
        else
            fInd = str2num(fn(1:end-4));
            fInd = sprintf('%04d', fInd);
        end
        copyfile([famImgs(fI).folder filesep famImgs(fI).name], [outPath filesep fInd '.jpg']);

    end
end

%% 

% path = 'G:\SUAnalysis\Fam_Task\P98CS\FamReScreenRecall_Session_1_20240803\match_by_grad_Nfam100_Nunfam120_1000';
% 
% ims = Utilities.readInFiles(path);
% ims = ims(1:100); % fam only
% 
% 
% for i = 1:length(ims)
%     fn = ims(i).name;
%     newN = 3000 + str2num(fn(1:end-4))
%     movefile([ims(i).folder filesep ims(i).name], [ims(i).folder filesep num2str(newN) '.jpg'])
% 
% end
% 
% 






