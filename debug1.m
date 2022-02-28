clear all
close all
clc

addpath('VirtualCamera/');


Xrange = 1 :  0.1 : 3
Yrange = 1 : 0.1 : 3;
camera_radius = 20;
image_noise_level = 0.01;

CAM =  VirtualCamera(Xrange, Yrange, camera_radius, image_noise_level);




% choose one imge
chc = 4;

chc = 2;

cam = CAM.GT_CameraMatrix{chc};
[K, R, C] = VirtualCamera.decomposeCameraMatrix(cam);

img_pts = CAM.ImgPoints(chc).Data;

img_norm = inv(K) * [img_pts; ones(1, length(img_pts))];
img_pts = img_norm([1,2], :) ./ img_norm(3, :);

template_pts = CAM.GT_Points3D;
template_pts = template_pts([1,2], :);


flag_training = false(1, size(img_pts, 2));
flag_training (1 : 7 : size(img_pts, 2)) = true;


fprintf(2, 'number of points for training: %f \n', sum(flag_training));
fprintf(2, 'number of points for validation: %f \n', sum(~flag_training));

hfig1 = figure(1);
hfig1.Position = [0, 0, 300, 300];
scatter(img_pts(1, flag_training), img_pts(2, flag_training), 'g*'); hold on;
scatter(img_pts(1, ~flag_training), img_pts(2, ~flag_training), 'k*'); hold on;
xlabel('x'); ylabel('y');

hfig2 = figure(2);
hfig2.Position = [400, 0, 300, 300];
scatter(template_pts(1, flag_training), template_pts(2, flag_training), 'go'); hold on;
scatter(template_pts(1, ~flag_training), template_pts(2, ~flag_training), 'ko'); hold on;
xlabel('x'); ylabel('y');


yvalues = img_pts(2, :);





JEstimate = ScanlineHomographyEstimation(img_pts(:, flag_training), template_pts(:, flag_training));
JEstimate.EstimateScanlineHomography();
scanlineHomographies = JEstimate.GetScanlineHomography(yvalues);
[poses, plane_homographies] = FundamentalHomographyEquation.GetScanlinePoses (yvalues, scanlineHomographies, 100,  'object');






GT_Pose = [R, -R*C]
pose_error = 0;
for ii = 1 : length(poses)    
    tmp = norm(poses{ii} - GT_Pose, 'fro');
    pose_error = pose_error + tmp*tmp;
end
pose_error = sqrt(pose_error/length(poses))




figure('Name', 'scanline homography mapping', 'Position', [0, 400, 600, 400]);
Jxpts_training = JEstimate.WarpRS2Template(img_pts(:, flag_training));
Jxpts_test = JEstimate.WarpRS2Template(img_pts(:, ~flag_training));
scatter(template_pts(1, flag_training), template_pts(2, flag_training), 'go'); hold on;
scatter(template_pts(1, ~flag_training), template_pts(2, ~flag_training), 'ko'); hold on;
scatter(Jxpts_training(1,:), Jxpts_training(2,:), 'g*'); 
scatter(Jxpts_test(1,:), Jxpts_test(2,:), 'k*');
hold off;
xlabel('x'); ylabel('y');
title(sprintf('Pose Error = %f', pose_error));




figure('Name', 'Scanline Homography', 'Position', [0, 1200, 1000, 200]); hold on;
Curve6 = zeros(6, length(scanlineHomographies));
for ii = 1 : length(scanlineHomographies)
    Curve6(:, ii) = reshape(scanlineHomographies{ii}, 6, 1);
end
[~, idx_order] = sort(yvalues);
k1 = plot(yvalues(idx_order),  Curve6(1, idx_order), '-', 'LineWidth', 2.5);
k2 = plot(yvalues(idx_order),  Curve6(2, idx_order), '-', 'LineWidth', 2.5);
k3 = plot(yvalues(idx_order),  Curve6(3, idx_order), '-', 'LineWidth', 2.5);
k4 = plot(yvalues(idx_order),  Curve6(4, idx_order), '-');
k5 = plot(yvalues(idx_order),  Curve6(5, idx_order), '-');
k6 = plot(yvalues(idx_order),  Curve6(6, idx_order), 'k-', 'LineWidth', 4);
hold off;
legend([k1, k2, k3, k4, k5, k6], {'k1', 'k2', 'k3', 'k4', 'k5', 'k6'});

