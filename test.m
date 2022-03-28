clear all
close all
clc


size_of_final_image = [NaN, 400];


anchor_idx = 1; 


fdir = '../CVPR_dataset/';

% this is the destination directory to contain figures for the paper
paper_figure_dir = '../PaperDraft/figure/';
if ~exist(paper_figure_dir, 'dir')
       mkdir(paper_figure_dir);
end

template_type = 'object';


data_options = [ 1 ];


params(1).param_paramterization_type = 'Polynomial';
params(1).param_polynomial_degree = [1, 1, 1, 1, 1];
params(1).param_num_control_points = [];

params(2).param_paramterization_type = 'Polynomial';
params(2).param_polynomial_degree = [1, 1, 1, 2, 2];
params(2).param_num_control_points = [];

params(3).param_paramterization_type = 'Polynomial';
params(3).param_polynomial_degree = [2, 2, 2, 2, 2];
params(3).param_num_control_points = [];

params(4).param_paramterization_type = 'Polynomial';
params(4).param_polynomial_degree = [2, 2, 2, 3, 3];
params(4).param_num_control_points = [];

params(5).param_paramterization_type = 'Polynomial';
params(5).param_polynomial_degree = [3, 3, 3, 3, 3];
params(5).param_num_control_points = [];

params(6).param_paramterization_type = 'BSpline';
params(6).param_polynomial_degree = [1, 1, 1, 1, 1];
params(6).param_num_control_points = [3, 3, 3, 3, 3];

params(7).param_paramterization_type = 'BSpline';
params(7).param_polynomial_degree = [1, 1, 1, 2, 2];
params(7).param_num_control_points = [3, 3, 3, 4, 4];

params(8).param_paramterization_type = 'BSpline';
params(8).param_polynomial_degree = [2, 2, 2, 2, 2];
params(8).param_num_control_points = [4, 4, 4, 4, 4];

params(9).param_paramterization_type = 'BSpline';
params(9).param_polynomial_degree = [2, 2, 2, 3, 3];
params(9).param_num_control_points = [4, 4, 4, 5, 5];

params(10).param_paramterization_type = 'BSpline';
params(10).param_polynomial_degree = [3, 3, 3, 3, 3];
params(10).param_num_control_points = [5, 5, 5, 5, 5];

params(11).param_paramterization_type = 'BSpline';
params(11).param_polynomial_degree = [3, 3, 3, 3, 3];
params(11).param_num_control_points = [5, 5, 5, 5, 5]+1;

params(12).param_paramterization_type = 'BSpline';
params(12).param_polynomial_degree = [3, 3, 3, 3, 3];
params(12).param_num_control_points = [5, 5, 5, 5, 5]+2;

params(13).param_paramterization_type = 'BSpline';
params(13).param_polynomial_degree = [3, 3, 3, 3, 3];
params(13).param_num_control_points = [5, 5, 5, 5, 5]+3;

params(14).param_paramterization_type = 'BSpline';
params(14).param_polynomial_degree = [3, 3, 3, 3, 3];
params(14).param_num_control_points = [5, 5, 5, 5, 5]+4;

params(15).param_paramterization_type = 'BSpline';
params(15).param_polynomial_degree = [3, 3, 3, 3, 3];
params(15).param_num_control_points = [5, 5, 5, 5, 5]+5;


params(16).param_paramterization_type = 'BSpline';
params(16).param_polynomial_degree = [3, 3, 3, 3, 3];
params(16).param_num_control_points = [5, 5, 5, 5, 5]+6;


params(17).param_paramterization_type = 'BSpline';
params(17).param_polynomial_degree = [3, 3, 3, 3, 3];
params(17).param_num_control_points = [5, 5, 5, 5, 5]+7;


params(18).param_paramterization_type = 'BSpline';
params(18).param_polynomial_degree = [3, 3, 3, 3, 3];
params(18).param_num_control_points = [5, 5, 5, 5, 5]+8;


for pcnt = 1 : length(params)

param_paramterization_type = params(pcnt).param_paramterization_type;    
param_polynomial_degree = params(pcnt).param_polynomial_degree;
param_num_control_points = params(pcnt).param_num_control_points;


for dtc = data_options

switch dtc

    case 1
        
        load('correspondence.mat');
        rollingshutter_name = 'RS';

    case 2

end

%%%%% -----------------------------

keypoints_template = keypoints.keypoints_template';
keypoints_rollingshutter = keypoints.keypoints_rolling_shutter';

calibration_rollingshutter = keypoints.intrinsics;
calibration_template = eye(3);

%%%%% -----------------------------------
image_rollingshutter = imread([rollingshutter_name, '.png']);


% solve Image rectification problem
RSPAPP = RollingShutterPlaneAbsolutePoseProblem;

RSPAPP.param_paramterization_type = param_paramterization_type;
RSPAPP.param_polynomial_degree = param_polynomial_degree;
RSPAPP.param_num_control_points = param_num_control_points;

[rectified_img, poses] = RSPAPP.SolveImageRectification(keypoints_rollingshutter, keypoints_template, calibration_rollingshutter, calibration_template, image_rollingshutter, template_type);
[rectified_landmarks,  lmk_poses] = RSPAPP.SolveLandmarkRectification (keypoints_rollingshutter, keypoints_template, calibration_rollingshutter, calibration_template, keypoints_rollingshutter, template_type);


anchor_idx = round(median(keypoints_rollingshutter(2, :)));

%%%%%  ground-truth rectification

% image rectification with ground-truth plane-homographies
gt_plane_homographies = keypoints.homography_per_rs_scanline;
gt_poses = keypoints.pose_per_rs_scaline;

height_t = size(keypoints.template,1);
step_t = height_t/length(gt_plane_homographies);


gt_rectified_img = ImageRectification.rectifyRSImageInverse(image_rollingshutter, gt_plane_homographies, anchor_idx);

imshow(gt_rectified_img); hold on;
grid on;
title('Complete rectified image with GT homographies');
hold off;
pause(0.1);

% landmark rectification with ground-truth
gt_lmk_poses = keypoints.global_shutter_pose_at_keypoints;

planarity_err = 0; mean_pose_err = 0;
backproj_pts = zeros(length(gt_lmk_poses),2);
gt_pts = zeros(length(gt_lmk_poses),2);


for ii = 1 : length(gt_lmk_poses)
    
    %%%%% OPTIONAL CODE (can be commented) %%%%%
    %%%%%%%% Code to test inverse projection of keypoints to template
    
% % % %     T =  gt_lmk_poses{ii};  % T = landmark poses
% % % %     T = [T; 0 0 0 1];
% % % %     T1 = T;
% % % %     R = T(1:3,1:3); % extracting rot. mat.
% % % %     t = T(1:3,4); % extracting translation vec.
% % % %     T = [R.' -(R.'*t); 0 0 0 1];  % T = T^{-1}, inverting transformation
% % % %     P = keypoints.keypoints_3D(ii,:).';  % P = 3D keypoints \in \mathbb{R}^3
% % % %     temp = P(2); % swapping X and Y axis
% % % %     P(2) = P(1);
% % % %     P(1) = temp;
% % % %     Pt = T*[P;1];  % Therefore: Pt = T^{-1} * P
% % % %     backproj_pts(ii,:) = Pt(1:2).';
% % % %     planarity_err = planarity_err + abs(Pt(3)); % any non-zero value in Z-axis is error, since the inverse transformed pointcloud should be perfecrtly planar and lying on XY-plane
% % % %     Pgt = keypoints.keypoints_template (ii,:); % Groundtruth template keypoints
% % % %     temp = Pgt(2); % Swapping axis for GT as well
% % % %     Pgt(2) = Pgt(1);
% % % %     Pgt(1) = temp;
% % % %     Pt = T*[P;1];
% % % %     gt_pts(ii,:) = Pgt(1:2).';
% % % %     mean_pose_err = mean_pose_err + ( (Pt(1) - Pgt(1))^2 + (Pt(2) - Pgt(2))^2); % computing mean pose error
% % % %     gt_lmk_poses{ii} = T1(1:3,:);
    
    %%%%%%%% Code to test inverse projection of keypoints to template
    %%%%% OPTIONAL CODE (can be commented) %%%%%


    gt_plane_homographies_lmks{ii} = (calibration_rollingshutter) * gt_lmk_poses{ii}(1:3, [1,2,4]);
end

%%%%% OPTIONAL CODE (can be commented) %%%%%
%%%%%%%% Code to test inverse projection of keypoints to template

% % % % scatter(backproj_pts(:,1),backproj_pts(:,2),'bo'); hold on; grid on;
% % % % scatter(gt_pts(:,1),gt_pts(:,2),'r*'); hold off;
% % % % legend('$\mathbf{T}^{-1}\mathbf{P}_{RS}$','$\mathbf{P}_{\tau}$','Interpreter','Latex');
% % % % title('$\mathbf{P}_{RS} = $RS pointcloud from kpts. in 3D; $\mathbf{P}_{\tau} = $ planar kpts. (GT) in object template','Interpreter','Latex');
% % % % drawnow();
% % % % pause(0.2);

%%%%%%%% Code to test inverse projection of keypoints to template
 %%%%% OPTIONAL CODE (can be commented) %%%%%

% the following code gives rectification of landmarks
H0 = gt_plane_homographies_lmks{anchor_idx};
for ii = 1 : length(keypoints_rollingshutter)
    pt = keypoints_rollingshutter(:, ii);
    H = gt_plane_homographies_lmks{ii};
    hgq = H0 * inv(H) * [pt; 1];
    gt_rectified_landmarks(:, ii) = hgq([1,2]) ./ hgq(3);
end

%%%%%

% swapT = [0 1 0 0; 1 0 0 0; 0 0 1 0; 0 0 0 1];
% for ii = 1 : length(gt_poses)
%     gt_poses{ii} = swapT*gt_poses{ii};
% end

[lmk_RMSE_rot, lmk_RMSE_pos] = internal_packages.TrajectoryError.RelativePoseError (lmk_poses, gt_lmk_poses);
[RMSE_rot, RMSE_pos] = internal_packages.TrajectoryError.RelativePoseError (poses, gt_poses);

translation_estimated = zeros(size(gt_poses,2), 3);
translation_gt = zeros(size(gt_poses,2), 3);

for ii = 1:size(gt_poses,2)
    translation_estimated(ii,:) = poses{ii}(1:3,4).';
    translation_gt(ii,:) = gt_poses{ii}(1:3,4).';
end

% close all;

figure;
subplot(3,1,1);
plot(translation_estimated(:,1),'k-','LineWidth',2); hold on;
plot(translation_gt(:,1),'r-','LineWidth',2); 
legend('$t_x$','$t_{x}^{GT}$','Interpreter','Latex');
hold off;
drawnow();

subplot(3,1,2);
plot(translation_estimated(:,2),'k-','LineWidth',2); hold on;
plot(translation_gt(:,2),'r-','LineWidth',2); 
legend('$t_y$','$t_{y}^{GT}$','Interpreter','Latex');
hold off;
drawnow();

subplot(3,1,3);
plot(translation_estimated(:,3),'k-','LineWidth',2); hold on;
plot(translation_gt(:,3),'r-','LineWidth',2); 
legend('$t_z$','$t_{z}^{GT}$','Interpreter','Latex');
hold off;
drawnow();

pause(1);

% close all;


pstr = [];
for k = 1 : length(RSPAPP.param_polynomial_degree)
    pstr = [pstr, num2str(RSPAPP.param_polynomial_degree(k))];
end
cstr = [];
for k = 1 : length(RSPAPP.param_num_control_points)
    cstr = [cstr, num2str(RSPAPP.param_num_control_points(k))];
end
if strcmp(RSPAPP.param_paramterization_type, 'Polynomial')
    append_info = ['polynomial_p', pstr];
end
if strcmp(RSPAPP.param_paramterization_type, 'BSpline')
    append_info = ['bspline_p', pstr, '_c', cstr];
end


if ~isnan(size_of_final_image(1))
scale_r = size_of_final_image(1)/size(image_rollingshutter, 1);
end
if ~isnan(size_of_final_image(2))
scale_r = size_of_final_image(2)/size(image_rollingshutter, 2);
end


normalized_keypoints_template = RSPAPP.output_KeyPoints.templatePointsNormalized;
template_Jxpts = RSPAPP.output_KeyPoints.Jx_rollingshutterPoints;

normalized_keypoints_template = RSPAPP.output_KeyPoints.templatePointsNormalized;
template_Jxpts = RSPAPP.output_KeyPoints.Jx_rollingshutterPoints;

figure('Name', 'Key Points Matching', 'Position', [0, 800, 500, 400]);
scatter(normalized_keypoints_template(1,:), normalized_keypoints_template(2,:), 'bo'); hold on;
scatter(template_Jxpts(1,:), template_Jxpts(2,:), 'r*'); hold off;
dvnorm = norm(template_Jxpts - normalized_keypoints_template, 'fro');
title(sprintf('RMSE = %f',  sqrt(dvnorm * dvnorm/ size(template_Jxpts, 2))));
xlabel('x'); ylabel('y');
pause(0.1);


%%% To be used with MATLAB 2019 or later
% % % % figure('Name', append_info, 'Position', [0, 0, 1900, 600]);
% % % % tfig = tiledlayout(1, 3, 'TileSpacing', 'tight');
% % % % ax1 = nexttile;
% % % % imshow(image_rollingshutter); pause(0.1); hold on;
% % % % scatter(keypoints_rollingshutter(1, :), keypoints_rollingshutter(2, :), '+', 'g'); pause(0.1);
% % % % hold off;
% % % % title(ax1, 'RS');
% % % % ax2 = nexttile;
% % % % imshow(rectified_img); hold on;
% % % % scatter(rectified_landmarks(1, :), rectified_landmarks(2, :), '+', 'g'); pause(0.1);
% % % % hold off;
% % % % title(ax2, 'Ours');
% % % % ax3 = nexttile;
% % % % imshow(gt_rectified_img); hold on;
% % % % scatter(gt_rectified_landmarks(1, :), gt_rectified_landmarks(2, :), '+', 'g'); pause(0.1);
% % % % hold off;
% % % % title(ax3, 'Ground-truth');
% % % % pause(1);


%%% To be used with MATLAB 2018 or earlier
figure('Name', 'RS', 'Position', [0, 100, 500, 400])
imshow(imresize(image_rollingshutter, size_of_final_image)); pause(0.1); hold on;
scatter(scale_r*keypoints_rollingshutter(1, :), scale_r*keypoints_rollingshutter(2, :), '+', 'g'); pause(0.1);
hold off;


figure('Name', 'Ours', 'Position', [800, 100, 500, 400])
imshow(imresize(rectified_img, size_of_final_image)); hold on;
scatter(scale_r*rectified_landmarks(1, :), scale_r*rectified_landmarks(2, :), '+', 'g'); pause(0.1);
hold off;


figure('Name', 'Ground-truth', 'Position', [800, 600, 500, 400]);
imshow(imresize(gt_rectified_img, size_of_final_image)); hold on;
scatter(scale_r*gt_rectified_landmarks(1, :), scale_r*gt_rectified_landmarks(2, :), '*', 'w'); pause(0.1);
title('Image rectified with GT homography');
hold off;


fprintf(2, '\nRelative Pose Error: \t scanline RMSE_pos = %f,    scanline RMSE_rot = %f\n', RMSE_pos, RMSE_rot);
fprintf(2, 'Relative Pose Error: \t keypoints RMSE_pos = %f,  keypoints RMSE_rot = %f\n\n', lmk_RMSE_pos, lmk_RMSE_rot);

pause(1);


end


end
