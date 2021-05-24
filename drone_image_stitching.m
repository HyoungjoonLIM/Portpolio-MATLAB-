%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%% Image Stitching %%%%%%%%%%%%%%%%%
%%%%%%%%%%%% 2018311004 Hyoung-joon Lim %%%%%%%%%%%
%%%%%%%%%%%%%%% In my lab-computer %%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all; clc

cd('E:\¼®»ç2ÇÐ±â\È«¹Ú»ç´Ô\image stiching')

image1 = imread('E:\¼®»ç2ÇÐ±â\È«¹Ú»ç´Ô\image stiching\1.jpg');
image2 = imread('E:\¼®»ç2ÇÐ±â\È«¹Ú»ç´Ô\image stiching\2.jpg');
image3 = imread('E:\¼®»ç2ÇÐ±â\È«¹Ú»ç´Ô\image stiching\3.jpg');
image4 = imread('E:\¼®»ç2ÇÐ±â\È«¹Ú»ç´Ô\image stiching\4.jpg');
totImages = 4;

IMAGE = [image1, image2, image3, image4];
IMAGE = im2double(rgb2gray(IMAGE));
ind = 7952*[1:totImages];

% Initialize features for I(1)
image1 = IMAGE(:,1:ind(1));
points = detectSURFFeatures(image1);
[features, points] = extractFeatures(image1, points);

% Initialize all the transforms to the identity matrix. Note that the
% projective transform is used here because the building images are fairly
% close to the camera. Had the scene been captured from a further distance,
% an affine transform would suffice.
tforms(totImages) = projective2d(eye(3));

% Initialize variable to hold image sizes.
imageSize = zeros(totImages,2);

% Iterate over remaining image pairs
for n = 2:totImages

    % Store points and features for I(n-1).
    pointsPrevious = points;
    featuresPrevious = features;

    % Read I(n).
    I = IMAGE(:,ind(n-1)+1:ind(n));

    % Save image size.
    imageSize(n,:) = size(I);

    % Detect and extract SURF features for I(n).
    points = detectSURFFeatures(I);
    [features, points] = extractFeatures(I, points);

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