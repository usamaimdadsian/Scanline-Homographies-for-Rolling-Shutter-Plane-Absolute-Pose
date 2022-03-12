function getDepthMap(frameLeftRect, frameRightRect)

frameLeftGray  = rgb2gray(frameLeftRect);
frameRightGray = rgb2gray(frameRightRect);

disparityRange = [0 64];
disparityMap = disparity(frameLeftGray, frameRightGray,'BlockSize',...
    15,'DisparityRange',disparityRange);
figure;
imshow(disparityMap, disparityRange);
title('Disparity Map');
colormap(gca,jet) 
colorbar

%% 3D reconsturction 
