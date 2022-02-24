clear all
close all
clc

% 
% 
% for i = 1 : 1000
% scanlineHomograph = rand(3,2);
% y = 10;
% 
% rota_trans = ScanlineHomograph.GetScanlinePose (y, scanlineHomograph);
% 
% 
% rota_trans = ScanlineHomograph.DisambiguatePose (rota_trans, [1, 0, 0]);
% 
% end

%VirtualCamera(Xrange, Yrange, camera_radius, image_noise_level);

Xrange = 1 : 20;
Yrange = 1 : 10;
camera_radius = 20;
image_noise_level = 0;
post_processing_noise = 0.00;

CAM =  VirtualCamera(Xrange, Yrange, camera_radius, image_noise_level);


% choose one imge
chc = 1;

cam = CAM.GT_CameraMatrix{chc};
[K, R, C] = VirtualCamera.decomposeCameraMatrix(cam);

img = CAM.ImgPoints(chc).Data;

img_norm = inv(K) * [img; ones(1, length(img))];
img = img_norm([1,2], :) ./ img_norm(3, :);


% figure(chc)
% scatter(img(1,:), img(2,:), 'r*'); hold on;
% xlabel('x'); ylabel('y');


nsize = length(Xrange);
% each line
for line = Yrange
imgeline = img(:, 1+(line-1)*nsize : nsize*line)
end

retinal_points = img.';
% img(:, [2,1]) = img;
template_points = CAM.GT_Points3D.';
template_points = template_points(:,1:2);
% retinal_scanline_ids = retinal_points(:,2);

retinal_noise = post_processing_noise*randn(size(retinal_points,1),2);
% retinal_points = retinal_points + retinal_noise;
% scatter(retinal_points(:,1),retinal_points(:,2),'b+');

rejection_identifier = -999999;
inlier_threshold = 0.004;
[horizon] = HorizontalLineExtraction(retinal_points, ones(size(retinal_points,1),1), 10, rejection_identifier, inlier_threshold);
tmp = 1 : length(horizon.scanline_ids);
ids = tmp(horizon.scanline_ids~=0);
D1 = retinal_points(horizon.scanline_ids~=rejection_identifier,1:2);
D2 = template_points(horizon.scanline_ids~=rejection_identifier,1:2);
ids = horizon.scanline_ids(horizon.scanline_ids~=rejection_identifier,1);

scatter(D1(:,1),D1(:,2),'g+'); hold on;
x = linspace(-0.5,0.5,20);
colmap = hsv(size(ids,1));
for ii = 1:size(ids,1)
    y = ones(20,1)*ids(ii);
    plot(x,y);
    hold on;
end
xlim([-0.5 0.5]);
ylim([-0.4 0.6]);
hold off;
% lines = [zeros(size(ids,1)) ids];
% scatter(D2(:,1),D2(:,2),'g+');

for ii = 1:size(ids,1)
    scatter(retinal_points(ii,1),retinal_points(ii,2),'bo');
    hold on;
    scatter(template_points(ii,1)+2,template_points(ii,2),'r+');
    plot([retinal_points(ii,1) template_points(ii,1)+2],[retinal_points(ii,2) template_points(ii,2)]);
end
hold off;

radialDistortion = [0.0 0.0 0.0];
intrinsics  = [K(1,1) 0 0; 0 K(2,2) 0; K(1,3), K(2,3) 1];
intrinsics = cameraParameters('IntrinsicMatrix',intrinsics,'RadialDistortion',radialDistortion);

data.retinal_points = D1;
data.template_points = D2;
data.retinal_scanline_ids = ids;
data.camera_parameters = intrinsics;
data.template_dimension = [0 0];
save('./../../Data/SynthMatlab_Fang_v1.7.mat', 'data');











