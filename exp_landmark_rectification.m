clear all
close all
clc



Xrange = 1 :  0.4 : 8;
Yrange = 1 : 2 : 10;
camera_radius = 20;
image_noise_level = 0.0;

rs_sigma_rot = 0.01;
rs_sigma_pos = 0.01;

%param_paramterization_type = 'Polynomial';
param_paramterization_type = 'BSpline';
param_polynomial_degree = [ 3, 3, 3, 3, 3];
param_num_control_points = [5, 5, 5, 5, 5];

CAM =  internal_packages.VirtualRollingShutterCamera (Xrange, Yrange, camera_radius, image_noise_level, rs_sigma_rot, rs_sigma_pos);



% choose one imge
chc = 1;

template_type = 'object';
% template_type = 'image';


rs_image = CAM.ImgPoints(chc);
rs_camera = CAM.GT_RollingShutterCameraMatrix(chc);

GT_poses = rs_camera.PerPointPoses;

img_pts = rs_image.RollingShutterAll;
GS_img_pts = rs_image.GlobalShutter;
calibration_rollingshutter = rs_camera.CalibrationMatrix;

flag_training = false(1, size(img_pts, 2));
flag_training (1 : 4 : size(img_pts, 2)) = true;
fprintf(1, 'number of points for training: %f \n', sum(flag_training));
fprintf(1, 'number of points for validation: %f \n', sum(~flag_training));


if strcmp(template_type, 'object')
    template_pts = CAM.GT_Points3D.All([1, 2], :);
    calibration_template = eye(3);
end
if strcmp(template_type, 'image')
    template_pts = CAM.ImgPoints(4).GlobalShutter;
    calibration_template = CAM.GT_RollingShutterCameraMatrix(4).CalibrationMatrix;
end






% solve landmark rectification problem
RSPAPP = RollingShutterPlaneAbsolutePoseProblem;
RSPAPP.param_paramterization_type = param_paramterization_type;
RSPAPP.param_polynomial_degree = param_polynomial_degree;
RSPAPP.param_num_control_points = param_num_control_points;

keypoints_rollingshutter = img_pts(:, flag_training);
keypoints_template = template_pts(:, flag_training);



landmarks_rollingshutter = img_pts;



[rectified_landmarks,  poses] = RSPAPP.SolveLandmarkRectification (keypoints_rollingshutter, keypoints_template, calibration_rollingshutter, calibration_template, landmarks_rollingshutter, template_type);




[APE_rot, APE_pos] = internal_packages.TrajectoryError.AbsoluteTrajectoryErrorByFirst(poses, GT_poses)
[RPE_rot, RPE_pos] = internal_packages.TrajectoryError.RelativePoseError(poses, GT_poses)
figure('Name', 'Trajectory Visualization', 'Position', [1400, 0, 600, 600]);
h1 = plotPoses (poses);
h2 = plotPoses(GT_poses);


Jxpts_landmarks = RSPAPP.output_RollingShutterLandmarks.Jx_rollingshutterLandmarks;
normalized_template_pts = inv(calibration_template) * [template_pts; ones(1, size(template_pts, 2));];
normalized_template_pts = normalized_template_pts([1,2], :) ./ normalized_template_pts(3, :);

figure('Name', 'scanline homography mapping', 'Position', [0, 400, 600, 400]);
hold on;
scatter(normalized_template_pts(1, flag_training), normalized_template_pts(2, flag_training), 'go');
scatter(normalized_template_pts(1, ~flag_training), normalized_template_pts(2, ~flag_training), 'ko');
scatter(Jxpts_landmarks(1, flag_training), Jxpts_landmarks(2, flag_training), 'g*');
scatter(Jxpts_landmarks(1, ~flag_training), Jxpts_landmarks(2, ~flag_training), 'k*');
hold off;
xlabel('x'); ylabel('y');





hfig1 = figure('Name', 'Rolling Shutter Landmarks', 'Position', [0, 0, 300, 300]);
hold on;
scatter(img_pts(1, flag_training), img_pts(2, flag_training), 'g*');
scatter(img_pts(1, ~flag_training), img_pts(2, ~flag_training), 'k*');
scatter(GS_img_pts(1, :),  GS_img_pts(2, :),  'rs');
hold off;
axis equal;
xlabel('x'); ylabel('y');


hfig2 = figure('Name', 'template', 'Position', [450, 0, 300, 300]);
hold on;
scatter(template_pts(1, flag_training), template_pts(2, flag_training), 'gs');
scatter(template_pts(1, ~flag_training), template_pts(2, ~flag_training), 'ks');
hold off;
axis equal;
xlabel('x'); ylabel('y');

hfig3 = figure('Name', 'Rectified Landmarks', 'Position', [800, 0, 300, 300]);
hold on;
scatter(rectified_landmarks(1, :), rectified_landmarks(2, :), '*');
scatter(GS_img_pts(1, :),  GS_img_pts(2, :),  'rs');
hold off;
axis equal;
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





function h = plotPoses (poseCell)
    pos_arr = [];
    for ii = 1 : length(poseCell)
        pos_arr = [pos_arr, poseCell{ii}(:, 4)];
    end
    hold on;
    h = plot3(pos_arr(1, :), pos_arr(2, :), pos_arr(3, :), '.-');
    hold off;
    view(3);
end