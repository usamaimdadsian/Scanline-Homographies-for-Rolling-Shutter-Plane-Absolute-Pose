clear all
close all
clc


size_of_final_image = [NaN, 400];


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
params(16).param_num_control_points = [5, 5, 5, 5, 5]+5;


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

keypoints_template = keypoints.keypoints_template(:, [2,1])';
keypoints_rollingshutter = keypoints.keypoints_rolling_shutter(:, [2,1])';

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



%[RMSE_rot, RMSE_pos] = internal_packages.TrajectoryError.RelativePoseError (lmk_poses, keypoints.global_shutter_pose_at_keypoints);


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




% [RMSE_rot, RMSE_pos] =
% internal_packages.TrajectoryError.RelativePoseError (poses, gt_poses);

close all;


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



figure('Name', append_info, 'Position', [0, 0, 1600, 600]);
tfig = tiledlayout(1, 2, 'TileSpacing', 'tight');
ax1 = nexttile;
imshow(image_rollingshutter); pause(0.1); hold on;
scatter(keypoints_rollingshutter(1, :), keypoints_rollingshutter(2, :), '.', 'r'); pause(0.1);
hold off;
ax2 = nexttile;
imshow(rectified_img); pause(0.1);



pause(1)

end


end
