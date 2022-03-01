clear all
close all
clc

addpath('VirtualCamera/');


Xrange = 1 :  0.1 : 3
Yrange = 1 : 1 : 5;
camera_radius = 20;
image_noise_level = 0.0;

rs_sigma_rot = 0.04;
rs_sigma_pos = 0.1;

param_paramterization_type = 'BSpline';
param_polynomial_degree = [ 3, 3, 3, 3, 3];
param_num_control_points = [5, 5, 5, 5, 5];

CAM =  internal_packages.VirtualRollingShutterCamera (Xrange, Yrange, camera_radius, image_noise_level, rs_sigma_rot, rs_sigma_pos);



% choose one imge
chc = 4;

chc = 2;


rs_image = CAM.ImgPoints(chc);
rs_camera = CAM.GT_RollingShutterCameraMatrix(chc);

K = rs_camera.CalibrationMatrix;
GT_poses = rs_camera.PerPointPoses;



img_pts = rs_image.All;
template_pts = CAM.GT_Points3D.All;


template_pts = template_pts([1,2], :);


flag_training = false(1, size(img_pts, 2));
flag_training (1 : 7 : size(img_pts, 2)) = true;


fprintf(2, 'number of points for training: %f \n', sum(flag_training));
fprintf(2, 'number of points for validation: %f \n', sum(~flag_training));

hfig1 = figure('Name', 'Rolling Shutter Landmarks', 'Position', [0, 0, 300, 300]);
scatter(img_pts(1, flag_training), img_pts(2, flag_training), 'g*'); hold on;
scatter(img_pts(1, ~flag_training), img_pts(2, ~flag_training), 'k*'); hold on;
xlabel('x'); ylabel('y');

hfig2 = figure('Name', 'template', 'Position', [400, 0, 300, 300]);
scatter(template_pts(1, flag_training), template_pts(2, flag_training), 'go'); hold on;
scatter(template_pts(1, ~flag_training), template_pts(2, ~flag_training), 'ko'); hold on;
xlabel('x'); ylabel('y');



% solve landmark rectification problem
RSPAPP = RollingShutterPlaneAbsolutePoseProblem;
RSPAPP.param_paramterization_type = param_paramterization_type;
RSPAPP.param_polynomial_degree = param_polynomial_degree;
RSPAPP.param_num_control_points = param_num_control_points;

keypoints_rollingshutter = img_pts(:, flag_training);
keypoints_template = template_pts(:, flag_training);

calibration_rollingshutter = K;
calibration_template = eye(3);

landmarks_rollingshutter = img_pts;
template_type = 'object';


[rectified_landmarks,  poses] = RSPAPP.SolveLandmarkRectification (keypoints_rollingshutter, keypoints_template, calibration_rollingshutter, calibration_template, landmarks_rollingshutter, template_type);





pose_error = 0;
for ii = 1 : length(poses)
    tmp = norm(poses{ii} - GT_poses{ii}, 'fro');
    pose_error = pose_error + tmp*tmp;
end
pose_error = sqrt(pose_error/length(poses))


Jxpts_training = RSPAPP.output_KeyPoints.Jx_rollingshutterPoints;
Jxpts_test = RSPAPP.output_RollingShutterLandmarks.Jx_rollingshutterLandmarks;

figure('Name', 'scanline homography mapping', 'Position', [0, 400, 600, 400]);
scatter(template_pts(1, flag_training), template_pts(2, flag_training), 'go'); hold on;
scatter(template_pts(1, ~flag_training), template_pts(2, ~flag_training), 'ko'); hold on;
scatter(Jxpts_training(1,:), Jxpts_training(2,:), 'g*'); 
scatter(Jxpts_test(1,:), Jxpts_test(2,:), 'k*');
hold off;
xlabel('x'); ylabel('y');
title(sprintf('Pose Error = %f', pose_error));






figure('Name', 'Rectified Landmarks', 'Position', [1000, 0, 700, 700]);
scatter(rectified_landmarks(1, :), rectified_landmarks(2, :), '*');
xlabel('x'); ylabel('y');




fontSize2 = 8;

scanlineHomographies = RSPAPP.output_RollingShutterLandmarks.scanlineHomographies;
yvalues = RSPAPP.output_RollingShutterLandmarks.yvaluesNormalized;

figure('Name', 'Scanline Homography', 'Position', [800, 1200, 1000, 300]);
tfig = tiledlayout(1, 2, 'TileSpacing','Loose');

ax1 = nexttile;
Curve6 = zeros(6, length(scanlineHomographies));
for ii = 1 : length(scanlineHomographies)
    Curve6(:, ii) = reshape(scanlineHomographies{ii}, 6, 1);
end
hold on;
k1 = plot(yvalues,  Curve6(1, :), '-', 'LineWidth', 2.5);
k2 = plot(yvalues,  Curve6(2, :), '-', 'LineWidth', 2.5);
k3 = plot(yvalues,  Curve6(3, :), '-', 'LineWidth', 2.5);
k4 = plot(yvalues,  Curve6(4, :), '-');
k5 = plot(yvalues,  Curve6(5, :), '-');
k6 = plot(yvalues,  Curve6(6, :), 'k-', 'LineWidth', 4);
hold off;
box on;
xlabel('normalized $y$ values', 'Interpreter','latex');
ylabel('$\mathbf{\Gamma}(y)$', 'Interpreter','latex');
title('un-normalized $\mathbf{J}(y)$', 'Interpreter','latex')

ax2 = nexttile;
Curve6 = zeros(6, length(scanlineHomographies));
for ii = 1 : length(scanlineHomographies)
    J = scanlineHomographies{ii};
    J = J ./ norm(J(:, 1), 'fro');
    Curve6(:, ii) = reshape(J, 6, 1);
end
hold on;
nk1 = plot(yvalues,  Curve6(1, :), '-', 'LineWidth', 2.5);
nk2 = plot(yvalues,  Curve6(2, :), '-', 'LineWidth', 2.5);
nk3 = plot(yvalues,  Curve6(3, :), '-', 'LineWidth', 2.5);
nk4 = plot(yvalues,  Curve6(4, :), '-');
nk5 = plot(yvalues,  Curve6(5, :), '-');
nk6 = plot(yvalues,  Curve6(6, :), 'k-', 'LineWidth', 4);
hold off;
box on;
xlabel('normalized $y$ values', 'Interpreter','latex');
ylabel('$\mathbf{\Gamma}(y)$', 'Interpreter','latex');
title('normalized $\mathbf{J}(y)$', 'Interpreter','latex')

lgd = legend([nk1, nk2, nk3, nk4, nk5, nk6], {'$\gamma_1$', '$\gamma_2$', '$\gamma_3$', '$\gamma_4$', '$\gamma_5$', '$\gamma_6$'}, ...
    'FontSize',fontSize2,'Interpreter','latex', 'Orientation', 'Horizontal', 'box', 'off');
lgd.Layout.Tile = 'South';


