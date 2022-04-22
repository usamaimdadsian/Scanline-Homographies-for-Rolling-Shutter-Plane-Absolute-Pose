clear all
close all
clc

addpath(genpath('./RollingShutterPlaneAbsolutePose'));

size_of_final_image = [NaN, 400];


anchor_idx = 1; 


% this is the destination directory to contain figures for the paper
paper_figure_dir = '../PaperDraft/figure/';
if ~exist(paper_figure_dir, 'dir')
       mkdir(paper_figure_dir);
end

template_type = 'object';


data_options = [ 1, 2,  3 ];


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



RMSE_rot_arr  = zeros(length(data_options), length(params));
RMSE_pos_arr = zeros(length(data_options), length(params));
RMSE_JxptFitness_arr = zeros(length(data_options), length(params));


exp_curve_names = {};



for pcnt = 1 : length(params)

param_paramterization_type = params(pcnt).param_paramterization_type;    
param_polynomial_degree = params(pcnt).param_polynomial_degree;
param_num_control_points = params(pcnt).param_num_control_points;


for dtc = data_options

switch dtc

    case 1

        fdir = 'data/RS29/';

        load([fdir, 'correspondence.mat']);
        rollingshutter_name = [fdir, 'RS'];
        dataName = 'ArbitraryMotionSim1';

        copyfile([fdir, 'op/R6P.png'],  [paper_figure_dir, dataName, '_R6P.png']);
    case 2

        fdir = 'data/RS30/';

        load([fdir, 'correspondence.mat']);
        rollingshutter_name = [fdir, 'RS'];
        dataName = 'ArbitraryMotionSim2';

        copyfile([fdir, 'op/R6P.png'],  [paper_figure_dir, dataName, '_R6P.png']);

    case 3

        fdir = 'data/RS31/';

        load([fdir, 'correspondence.mat']);
        rollingshutter_name = [fdir, 'RS'];
        dataName = 'ArbitraryMotionSim3';

        copyfile([fdir, 'op/R6P.png'],  [paper_figure_dir, dataName, '_R6P.png']);

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


gt_rectified_img = ImageRectification.rectifyRSImage (image_rollingshutter, gt_plane_homographies, anchor_idx);

% imshow(gt_rectified_img); hold on;
% grid on;
% title('Complete rectified image with GT homographies');
% hold off;
% pause(0.1);

% landmark rectification with ground-truth
gt_lmk_poses = keypoints.global_shutter_pose_at_keypoints;

% planarity_err = 0; mean_pose_err = 0;
% backproj_pts = zeros(length(gt_lmk_poses),2);
% gt_pts = zeros(length(gt_lmk_poses),2);


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

% [lmk_RMSE_rot, lmk_RMSE_pos] = internal_packages.TrajectoryError.AbsoluteTrajectoryErrorByFirst(lmk_poses, gt_lmk_poses);
% [RMSE_rot, RMSE_pos] = internal_packages.TrajectoryError.AbsoluteTrajectoryErrorByFirst(poses, gt_poses);

[lmk_RMSE_rot, lmk_RMSE_pos] = internal_packages.TrajectoryError.RelativePoseError (lmk_poses, gt_lmk_poses);
[RMSE_rot, RMSE_pos] = internal_packages.TrajectoryError.RelativePoseError (poses, gt_poses);

RMSE_rot_arr(dtc, pcnt) = RMSE_rot;
RMSE_pos_arr(dtc, pcnt) = RMSE_pos;



translation_estimated = zeros(size(gt_poses,2), 3);
translation_gt = zeros(size(gt_poses,2), 3);

for ii = 1:size(gt_poses,2)
    translation_estimated(ii,:) = poses{ii}(1:3,4).';
    translation_gt(ii,:) = gt_poses{ii}(1:3,4).';
end

% close all;

% figure;
% subplot(3,1,1);
% plot(translation_estimated(:,1),'k-','LineWidth',2); hold on;
% plot(translation_gt(:,1),'r-','LineWidth',2); 
% legend('$t_x$','$t_{x}^{GT}$','Interpreter','Latex');
% hold off;
% drawnow();
% 
% subplot(3,1,2);
% plot(translation_estimated(:,2),'k-','LineWidth',2); hold on;
% plot(translation_gt(:,2),'r-','LineWidth',2); 
% legend('$t_y$','$t_{y}^{GT}$','Interpreter','Latex');
% hold off;
% drawnow();
% 
% subplot(3,1,3);
% plot(translation_estimated(:,3),'k-','LineWidth',2); hold on;
% plot(translation_gt(:,3),'r-','LineWidth',2); 
% legend('$t_z$','$t_{z}^{GT}$','Interpreter','Latex');
% hold off;
% drawnow();

% pause(1);

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
    append_name_info = ['p', pstr];
end
if strcmp(RSPAPP.param_paramterization_type, 'BSpline')
    append_info = ['bspline_p', pstr, '_c', cstr];
    append_name_info = ['p', pstr, '\_c', cstr];
end
exp_curve_names = [exp_curve_names,  append_name_info];



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


RMSE_JxptFitness_arr(dtc, pcnt) = sqrt(dvnorm * dvnorm/ size(template_Jxpts, 2));

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
rsfig = figure('Name', 'RS', 'Position', [0, 100, 500, 400]);
imshow(imresize(image_rollingshutter, size_of_final_image)); pause(0.1); hold on;
scatter(scale_r*keypoints_rollingshutter(1, :), scale_r*keypoints_rollingshutter(2, :), '+', 'g'); pause(0.1);
hold off;
exportgraphics(rsfig, [paper_figure_dir, dataName, '_rs_img.pdf'], 'ContentType', 'vector');
% imwrite(image_rollingshutter, [paper_figure_dir, dataName, '_rs_img.png']);


rectfig = figure('Name', ['Ours: ', append_info], 'Position', [800, 100, 500, 400]);
imshow(imresize(rectified_img, size_of_final_image)); hold on;
scatter(scale_r*rectified_landmarks(1, :), scale_r*rectified_landmarks(2, :), '+', 'g'); pause(0.1);
hold off;
exportgraphics(rectfig, [paper_figure_dir, dataName, '_rect_img_', append_info, '.pdf'], 'ContentType', 'vector');
% imwrite(rectified_img, [paper_figure_dir, dataName, '_rect_img_', append_info, '.png']);


gtfig = figure('Name', 'Ground-truth Rectification', 'Position', [800, 600, 500, 400]);
imshow(imresize(gt_rectified_img, size_of_final_image)); hold on;
scatter(scale_r*gt_rectified_landmarks(1, :), scale_r*gt_rectified_landmarks(2, :), '+', 'g'); pause(0.1);
hold off;
exportgraphics(gtfig, [paper_figure_dir, dataName, '_ground_truth_rectification.pdf'], 'ContentType', 'vector');
%imwrite(gt_rectified_img, [paper_figure_dir, dataName, '_ground_truth_rectification.png']);


fprintf(1, '\nRelative Pose Error: \t scanline RMSE_pos = %f,    scanline RMSE_rot = %f\n', RMSE_pos, RMSE_rot);
fprintf(1, 'Relative Pose Error: \t keypoints RMSE_pos = %f,  keypoints RMSE_rot = %f\n\n', lmk_RMSE_pos, lmk_RMSE_rot);

pause(1);

close all;

end


end


exp_curve_names = exp_curve_names(1 : length(data_options): end);

RMSE_pos_arr

RMSE_rot_arr


RMSE_JxptFitness_arr




%%

if (1)

    close all;


    datasetNames = {'RS1', 'RS2', 'RS3'};

    fontSize1 = 8;

    fontSize2 = 7;

    fontSize3 = 8;

    figure('Name', 'Image rectification', 'Position', [0, 800, 800, 250]);
    tfig = tiledlayout(1, 3, 'TileSpacing', 'compact');

    linespec = {'.-', '.-', '.-', '.-', '.-',     'o--', 'o--', 'o--', 'o--', 'o--'};
    linewith = {1.5, 1.5, 1.5, 1.5, 1.5,     2, 2, 2, 2, 2};
    linecolor = {'#A2142F', '#77AC30',  '#0072BD', '#4DBEEE', '#EDB120'};

    colors = colorcube(size(RMSE_rot_arr,2));

    ax1 = nexttile;
    hold on;
    b1 = bar(RMSE_rot_arr, 'FaceColor', 'flat');
    hold off; grid off; box on;
    for k = 1:size(RMSE_rot_arr,2)
    b1(k).CData = colors(k, :);
    end
    ylabel(ax1, 'Rotation RMSE', 'FontSize',fontSize1);
    ax1.FontSize = fontSize3;
    xticks(1 : length(data_options));
    xticklabels(datasetNames);


    ax2 = nexttile;
    hold on;
    b2 = bar(RMSE_pos_arr, 'FaceColor', 'flat');
    hold off; grid off; box on;
    for k = 1:size(RMSE_pos_arr,2)
        b2(k).CData = colors(k, :);
    end
    ylabel(ax2, 'Translation RMSE', 'FontSize',fontSize1);
    ax2.FontSize = fontSize3;
    xticks(1 : length(data_options));
    xticklabels(datasetNames);


    ax3 = nexttile;
    hold on;
    b2 = bar(RMSE_JxptFitness_arr, 'FaceColor', 'flat');
    hold off; grid off; box on;
    for k = 1:size(RMSE_JxptFitness_arr,2)
        b2(k).CData = colors(k, :);
    end
    %ylabel(ax3, '$\mathbf{J}(y)$ Fitness', 'Interpreter','latex', 'FontSize',fontSize1);
    ylabel(ax3, 'Fitness of J(y)', 'FontSize',fontSize1);
    ax3.FontSize = fontSize3;
    xticks(1 : length(data_options));
    xticklabels(datasetNames);

    lgd = legend(exp_curve_names, 'FontSize',fontSize2,  'Orientation','Vertical', 'box', 'off');

    lgd.Layout.Tile = 'east';


    pause(0.1);


    exportgraphics(tfig, [paper_figure_dir, 'ScanlinePoseError_ArbitraryMotionSimulation', '.pdf'], 'ContentType', 'Vector');

    return;

end



if (1)
    close all;

    linespec = {'.-', '.-', '.-', '.-', '.-',     'o--', 'o--', 'o--', 'o--', 'o--'};
    linewith = {1.2, 1.2, 1.2, 1.2, 1.2,     2, 2, 2, 2, 2};
    linecolor = {'#A2142F', '#77AC30',  '#0072BD', '#4DBEEE', '#EDB120'};

    Trans = [];
    gtTrans = [];
    for ii = 1 : length(poses)
        Trans = [Trans,  poses{ii}(1:3, 4)];
        gtTrans = [gtTrans, gt_poses{ii}(1:3, 4)];
    end

    figure('Name', 'Translation Error', 'Position', [0, 0, 500, 200]);
    tfig = tiledlayout(1, 3, 'TileSpacing', 'compact');

    ax1 = nexttile;
    hold on;
    plot(Trans(1, :), linespec{1}, 'LineWidth', linewith{1});
    plot(gtTrans(1, :), linespec{2}, 'LineWidth', linewith{2});
    hold off; box on;
    ylabel('X axis', 'Fontsize', 8);
    xlabel('scanline indices', 'Fontsize', 8)

    ax2 = nexttile;
    hold on;
    plot(Trans(2, :), linespec{1}, 'LineWidth', linewith{1});
    plot(gtTrans(2, :), linespec{2}, 'LineWidth', linewith{2});
    hold off; box on;
    ylabel('Y axis', 'Fontsize', 8);
    xlabel('scanline indices', 'Fontsize', 8)

    ax3 = nexttile;
    hold on;
    plot(Trans(3, :), linespec{1}, 'LineWidth', linewith{1});
    plot(gtTrans(3, :), linespec{2}, 'LineWidth', linewith{2});
    hold off; box on;
    ylabel('Z axis', 'Fontsize', 8);
    xlabel('scanline indices', 'Fontsize', 8)


    lgd = legend({'Estimated', 'Ground-truth'}, 'FontSize', 8,'Interpreter','latex', 'box', 'off', 'Orientation','horizontal');
    lgd.Layout.Tile = 'south';   

    exportgraphics(tfig, [paper_figure_dir, 'ScanlinePoseError_ArbitraryMotionSimulation_XYZ', '.pdf'], 'ContentType', 'Vector');

end





