clear all
close all
clc

addpath(genpath('./RollingShutterPlaneAbsolutePose'));

Xrange = 1 :  0.2 : 6;
Yrange = 1 : 2 : 10;
camera_radius = 20;
image_noise_level = 0.0;

rs_sigma_rot = 0.01;
rs_sigma_pos = 0.1;

%param_paramterization_type = 'Polynomial';
param_paramterization_type = 'BSpline';
param_polynomial_degree = [ 2, 2, 2, 3, 3];
param_num_control_points = [4, 4, 4, 5, 5];

CAM =  internal_packages.VirtualRollingShutterCamera (Xrange, Yrange, camera_radius, image_noise_level, rs_sigma_rot, rs_sigma_pos);



% choose one imge
chc = 1;

template_type = 'object';
%template_type = 'image';


rs_image = CAM.ImgPoints(chc);
rs_camera = CAM.GT_RollingShutterCameraMatrix(chc);

GT_poses = rs_camera.PerPointPoses;

img_pts = rs_image.RollingShutterAll;
GS_img_pts = rs_image.GlobalShutter;
calibration_rollingshutter = rs_camera.CalibrationMatrix;

flag_training = false(1, size(img_pts, 2));
flag_training (1 : 5 : size(img_pts, 2)) = true;
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



figure('Name', 'Scene visualization')
hold on;
scatter3(CAM.GT_Points3D.All(1, :), CAM.GT_Points3D.All(2, :), CAM.GT_Points3D.All(3, :), 'ks')
camera_pose = rs_camera.PerScanlinePoses{CAM.global_shutter_anchor};
% cam = PlotCameraByPose (camera_pose);
hold off;
view(3);
% cam.AxesVisible = true;



% solve landmark rectification problem
RSPAPP = RollingShutterPlaneAbsolutePoseProblem;
RSPAPP.param_paramterization_type = param_paramterization_type;
RSPAPP.param_polynomial_degree = param_polynomial_degree;
RSPAPP.param_num_control_points = param_num_control_points;

keypoints_rollingshutter = img_pts(:, flag_training);
keypoints_template = template_pts(:, flag_training);



landmarks_rollingshutter = img_pts;



[rectified_landmarks,  poses] = RSPAPP.SolveLandmarkRectification (keypoints_rollingshutter, keypoints_template, calibration_rollingshutter, calibration_template, landmarks_rollingshutter, template_type);


[rectified_landmarks_RSPnP_Velocity,  poses_RSPnP_Velocity] = PlanarRSPnP.SolveLandmarkRectification (keypoints_rollingshutter, keypoints_template, calibration_rollingshutter, calibration_template, landmarks_rollingshutter, template_type);



if (0)
    
    
    [dummy_image, cX, cY] = PrepDummyImageForYizhen(keypoints_rollingshutter, keypoints_template, 500); % blank images for Yizhen's method
    
    calibration_rollingshutter(1,3) = cX; calibration_rollingshutter(2,3) = cY; % updating principal points in 'K'
    calibration_template = calibration_rollingshutter; % template and RS image needs to necessarily have same intrinsics for Yizhen
    
    yRS = BenchmarkYizhen; % creating the object
    keypoints_rollingshutter = keypoints_rollingshutter + [cX;cY]; % pushing keypoints away from center
    keypoints_template = keypoints_template + [cX;cY]; % pushing keypoints away from center
    landmarks_rollingshutter = landmarks_rollingshutter + [cX;cY]; % pushing landmarks away from center
    
    %%%%%% Plotting input %%%%%%
    scatter(keypoints_rollingshutter(1,:),keypoints_rollingshutter(2,:),'r*'); hold on;
    scatter(keypoints_template(1,:),keypoints_template(2,:),'b*');
    legend('RS keypoints','Template Keypoints','FontSize',14);
    title('RS and template keypoints supplied to Yizhens method','FontSize',12);
    grid on; hold off; drawnow(); pause(1);
    %%%%%% Plotting input %%%%%%
    
    % Yizhen's method
    [rectified_image_Yizhen, landmark_poses_Yizhen, landmark_homography_Yizhen, status] = yRS.SolveImageRectification(keypoints_rollingshutter, keypoints_template, calibration_rollingshutter, dummy_image, dummy_image, landmarks_rollingshutter); % firing Yizhen's method
    
    % landmark poses rectified (obtained by applying landmark homographies (???))
    landmark_poses_Yizhen = landmark_poses_Yizhen.';
    landmark_poses_Yizhen = landmark_poses_Yizhen./landmark_poses_Yizhen(3,:); % projecting
    landmark_poses_Yizhen(3,:) = []; % discarding 3rd column
    landmark_poses_Yizhen = landmark_poses_Yizhen - [cX;cY]; % pulling back the landmark points towards center, offsetting the displacement by principal points
    
    %%%%% Plotting the output %%%%
    scatter(rectified_landmarks(1,:), rectified_landmarks(2,:),'r+'); hold on
    scatter(landmark_poses_Yizhen(1,:), landmark_poses_Yizhen(2,:),'b*');
    legend('Rectified landmarks: OUR method','Rectified landmarks: YIZHEN method','FontSize',12);
    grid on; hold off;
    %%%%% Plotting the output %%%%
    
    
end


[APE_rot, APE_pos] = internal_packages.TrajectoryError.AbsoluteTrajectoryErrorByFirst(poses, GT_poses)
[RPE_rot, RPE_pos] = internal_packages.TrajectoryError.RelativePoseError(poses, GT_poses)
% figure('Name', 'Trajectory Visualization', 'Position', [1400, 0, 600, 600]);




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
% tfig = tiledlayout(1, 2, 'TileSpacing','Loose');

% ax1 = nexttile;
subplot(1,2,1);
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

% ax2 = nexttile;
subplot(1,2,2);
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
% lgd.Layout.Tile = 'South';

pause(2);







function [image, cX, cY] = PrepDummyImageForYizhen(keypoints_rollingshutter, keypoints_template, buffer_pix)
    keypoints_fused = [keypoints_rollingshutter keypoints_template];
    
    maxX = ceil(max(keypoints_fused(1,:))) + buffer_pix;
    maxY = ceil(max(keypoints_fused(2,:))) + buffer_pix;
    
    cX = maxX/2;
    cY = maxY/2;
    
    image = ones(maxX, maxY,3);
    
end


function cam = PlotCameraByPose (pose)
        rotationMatrix = pose(1:3, 1:3);
        translationVector = pose(1:3, 4)';
        orientation = rotationMatrix';
        location = -translationVector * orientation;        
        absPose = rigid3d(orientation,location);
        cam = plotCamera('AbsolutePose',absPose);
end


