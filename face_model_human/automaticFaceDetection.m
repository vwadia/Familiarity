%  Using matlab's built in computer vision toolbox to detect face landmarks
%  maybe this will work better than the python implementation

run G:\SUAnalysis\setDiskPaths

% read into photos
imPath = [famPath filesep 'LiangsFamPhotos'];
ims = Utilities.readInFiles(imPath, 'jpg');
ims = [ims(1); ims(2)];

for i = 1:length(ims)
    
    img = imread([ims(i).folder filesep ims(i).name]);

    % this only has 5 landmarks
    [bboxes, scores, landmarks] = mtcnn.detectFaces(img);
    imshow(img)
    
    faceDetector = vision.CascadeObjectDetector();
    bboxes = step(faceDetector, img);
    figure
    imshow(img)
    
end

%% 
data = importdata('shape_predictor_68_face_landmarks.dat');