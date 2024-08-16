function pts = get_landmark92(file80, toSave)
% Bastardized version ot Liangs interpolation code 
% takes in a struct with 80 points, interpolates to get the last 12
% the returns the new set of 92 points
% 
% INPUTS:
%     1. Struct pointing to file
%     
% OUTPATH:
%     1. 92x2 pts double
%
% vwadia May2023

if nargin == 1, toSave = false; end

[path, filename, ext] = fileparts([file80.folder filesep file80.name]);

file92 = fullfile(path, [filename ext]);

load([file80.folder filesep file80.name]);
pts = double(pts);

% don't add extra points
if length(pts) == 80
    pt_ind = [1 79 74 76 78 80 17];
    pti = pts(pt_ind, :); % pt to interpolate
    
    xq = setdiff(1:1/3: 7, 1:7)';
    pt_x = interp1(pti(:,1), xq, 'pchip');
    pt_y = interp1(pti(:,2), xq, 'pchip');
    
    pts = [pts; [pt_x pt_y]];
    
    if toSave
        save(file92, 'pts');
    end
end
   
end
