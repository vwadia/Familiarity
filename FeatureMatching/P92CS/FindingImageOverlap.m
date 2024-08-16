% Script to make sure I don't have to re-mark all the points for images I
% am re-using
% vwadia 2024


% P98 images
P98ImPath = 'C:\Users\wadiav\Documents\PYTHON\FaceFamiliarity\BackgroundRemoval_SegAny\P98CS\BigSet_Processed';

% P92 images
P92ImPath = 'C:\Users\wadiav\Documents\PYTHON\FaceFamiliarity\BackgroundRemoval_SegAny\P92CS\FamFaces_Pt_P92CS_Processed';


p98Ims = Utilities.readInFiles(P98ImPath);
p92Ims = Utilities.readInFiles(P92ImPath);

smvals = nan(length(p98Ims), length(p92Ims), 3);

for i = 1:length(p98Ims)
    targ = imread([p98Ims(i).folder filesep p98Ims(i).name]);
    % tic
    for j = 1:length(p92Ims)

        check = imread([p92Ims(j).folder filesep p92Ims(j).name]);
        smvals(i, j, :) = multissim(targ, check);
        if mean(smvals(i, j, :)) == 1
            break
        end
                

    end
    % toc
    fprintf('Finished for P98 image number %d \n', i)
end

smmean = mean(smvals, 3);
smmean = smmean(1:end-1, :);

matchNums = nan(length(p98Ims)-1, 1);
% now find indices
for i = 1:length(p98Ims)-1
    if ~isempty(find(smmean(i, :) == 1))
        matchNums(i, 1) = find(smmean(i, :) == 1);
    end
end

% load in familiar indices
load([P98ImPath filesep 'ratings_1To3_Equals_UnfamToFam.mat']) % creates ids

IdsToSearch = cellfun(@(x) x > 2, ids(:, 1));

nonanNums = matchNums(~isnan(matchNums) & IdsToSearch);

% move Images
destFold = 'C:\Users\wadiav\Documents\PYTHON\FaceFamiliarity\BackgroundRemoval_SegAny\P98CS\FamFaces_P98CS_Processed';
ctr = 1;
for i = nonanNums'   
    fnum = sprintf('%04d', ctr);
    copyfile([p92Ims(i).folder filesep p92Ims(i).name], [destFold filesep fnum '.jpg'])
    ctr = ctr+1;
end


%% doing this by hand...
% copy these ims to folder and these pts to folder

% p92idsHand = [1 2 10:13 16 18 19 21 22 26 29 31:36 39 40 44 45 47 50 51 52 53 55 56 59:61 65 66 69 71 72 74:77 80 86 91];
% 
% p87idsHand = [3 59];
% 
% p86idsHand = [10 20 46 56 57 66 71 79 84:86]; % copy these ims to folder and these pts to folder
basePath = 'C:\Users\wadiav\Documents\PYTHON\FaceFamiliarity\BackgroundRemoval_SegAny';

sessID{1} = [basePath filesep 'P92CS' filesep 'FamFaces_Pt_P92CS_Processed'];
sessID{2} = [basePath filesep 'P87CS' filesep 'FamFaces_Pt_P87CS_2_Processed'];
sessID{3} = [basePath filesep 'P86CS' filesep 'FamFaces_Pt_P86CS_Processed'];

targSess = [basePath filesep 'P98CS' filesep 'toKeepFaces'];
targIms = Utilities.readInFiles(targSess);
sessIdx = [];

% for each session
for ss = 1:length(sessID)
    Ims = Utilities.readInFiles(sessID{ss});
    idx = ones(length(Ims), 1)+(ss-1);
    if ss ==1
        compIms = Ims;
    else
        compIms = cat(1, compIms,  Ims);
    end
    sessIdx = cat(1, sessIdx, idx);
end

compImages = cell(length(compIms), 1);
for c = 1:length(compIms)
    compImages{c} = imread([compIms(c).folder filesep compIms(c).name]);
end

destFoldPts = [basePath filesep 'P98CS' filesep 'markedPts_Pt_P98CS'];
if ~exist(destFoldPts)
    mkdir(destFoldPts)
end

destFoldIms = [basePath filesep 'P98CS' filesep 'FamFaces_P98CS_Processed'];
if ~exist(destFoldIms)
    mkdir(destFoldIms)
end


for t = 1:length(targIms)
    targ = imread([targIms(t).folder filesep targIms(t).name]);
    fnum = sprintf('%04d', t);
    for c = 1:length(compImages)
        smval = multissim(targ, compImages{c}); 
        if mean(smval) == 1
            % move image
            copyfile([compIms(c).folder filesep compIms(c).name], [destFoldIms filesep fnum '.jpg'])

            % move appropriate marked pts
            if sessIdx(c) == 1
                dfoldMP = [basePath filesep 'P92CS' filesep 'markedPts_Pt_P92CS'];
                cpts = c;
            elseif sessIdx(c) == 2
                dfoldMP = [basePath filesep 'P87CS' filesep 'markedPts_Pt_P87CS_2'];
                cpts = c - sum(sessIdx == 1);
            elseif sessIdx(c) == 3
                dfoldMP = [basePath filesep 'P86CS' filesep 'markedPts_Pt_P86CS'];
                cpts = c - (sum(sessIdx == 1) + sum(sessIdx == 2));
            end
            mpts = Utilities.readInFiles(dfoldMP);
            assert(isequal(str2num(compIms(c).name(1:end-4)), str2num(mpts(cpts).name(1:end-4))))

            copyfile([dfoldMP filesep mpts(cpts).name(1:end-4) '.mat'], [destFoldPts filesep fnum '.mat'])
        end
    end
    fprintf('Finished for image %d/%d \n', t, length(targIms))
end

%% 

basePath = 'C:\Users\wadiav\Documents\PYTHON\FaceFamiliarity\BackgroundRemoval_SegAny';
% imPath = [basePath filesep 'P98CS' filesep 'ExtraFaces_Raw_oldNames'];
% imPath = [basePath filesep 'P98CS' filesep 'FamFaces_P98CS_Processed'];
imPath = [basePath filesep 'P98CS' filesep 'FamFaces_P98CS_Raw'];
% imPath = [basePath filesep 'P98CS' filesep 'markedPts_Pt_P98CS'];

ims = Utilities.readInFiles(imPath);
% dFol = [basePath filesep 'P98CS' filesep 'FamFaces_P98CS_Raw'];
% dFol = [basePath filesep 'P98CS' filesep 'ExtraFaces_Raw'];
% dFol =  [basePath filesep 'P98CS' filesep 'NewIds_Im'];
dFol =  [basePath filesep 'P98CS' filesep 'NewIds_Im_Raw'];
% dFol =  [basePath filesep 'P98CS' filesep 'NewIds_Pts'];

for i = 1:length(ims)
     fnum = sprintf('%04d', i);
    copyfile([ims(i).folder filesep ims(i).name], [dFol filesep fnum '.jpg'])
end
