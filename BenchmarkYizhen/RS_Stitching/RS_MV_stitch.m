% Real RS images stitching  
% Author: Yizhen LAO
% Mutiple views example
% version: 1.0  
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc 
close all
clear all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Step 1 - Load Images
% As seen below, the images are relatively unaffected by any lens distortion so camera calibration was not required. 
% However, if lens distortion is present, the camera should be calibrated and the images undistorted prior to creating the panorama. 
% You can use the cameraCalibrator App to calibrate a camera if needed.
buildingDir = fullfile(toolboxdir('vision'), 'visiondata', 'building');
buildingScene = imageDatastore(buildingDir);

% Display images to be stitched
montage(buildingScene.Files)

%% Step 2 - Register Image Pairs
% To create the panorama, start by registering successive image pairs using the following procedure:
% 1. Detect and match features between  and .
% 2. Estimate the geometric transformation, , that maps  to .
% 3. Compute the transformation that maps  into the panorama image as .

% Read the first image from the image set.
I = readimage(buildingScene, 1);

% Initialize features for I(1)
grayImage = rgb2gray(I);
points = detectSURFFeatures(grayImage);
[features, points] = extractFeatures(grayImage, points);

% Initialize all the transforms to the identity matrix. Note that the
% projective transform is used here because the building images are fairly
% close to the camera. Had the scene been captured from a further distance,
% an affine transform would suffice.
numImages = numel(buildingScene.Files);
tforms(numImages) = projective2d(eye(3));

% Initialize variable to hold image sizes.
imageSize = zeros(numImages,2);

% Iterate over remaining image pairs
for n = 2:numImages

    % Store points and features for I(n-1).
    pointsPrevious = points;
    featuresPrevious = features;

    % Read I(n).
    I = readimage(buildingScene, n);

    % Convert image to grayscale.
    grayImage = rgb2gray(I);

    % Save image size.
    imageSize(n,:) = size(grayImage);

    % Detect and extract SURF features for I(n).
    points = detectSURFFeatures(grayImage);
    [features, points] = extractFeatures(grayImage, points);

    % Find correspondences between I(n) and I(n-1).
    indexPairs = matchFeatures(features, featuresPrevious, 'Unique', true);

    matchedPoints = points(indexPairs(:,1), :);
    matchedPointsPrev = pointsPrevious(indexPairs(:,2), :);

    % Estimate the transformation between I(n) and I(n-1).
    tforms(n) = estimateGeometricTransform(matchedPoints, matchedPointsPrev,...
        'projective', 'Confidence', 99.9, 'MaxNumTrials', 2000);

    % Compute T(n) * T(n-1) * ... * T(1)
    tforms(n).T = tforms(n).T * tforms(n-1).T;
end