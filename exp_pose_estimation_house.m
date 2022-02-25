clear all
close all
clc


size_of_final_image = [NaN, 400];

fdir = '../CVPR_dataset/HouseSimulation/';

subdir = {'house_rot1_B0',  'house_rot1_B20',  'house_rot1_B40',  'house_rot2_B40',  'house_trans1_B40',  'house_trans_rot1_B40'};
frame_names = {'frame01', 'frame02', 'frame03', 'frame04', 'frame05', 'frame06', 'frame07', 'frame08', 'frame09', 'frame10', 'frame11', 'frame12'};
calibration_matrix = readmatrix([fdir, 'CalibrationMatrix', '.txt']);

data_options =[6];
data_choice = data_options(1);
data_dir = [fdir, subdir{data_choice}];

%template_image = 'gt_start_frame01_template.png';
%template_image = 'gt_start_frame09_template.png'
%gs_template_image = imread(template_image);




% this is the destination directory to contain figures for the paper
paper_figure_dir = '../PaperDraft/figure/';
if ~exist(paper_figure_dir, 'dir')
       mkdir(paper_figure_dir);
end



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



% pcnt = [1, 2, 3]
RMSE_rot_arr = [];
RMSE_pos_arr = [];

exp_curve_names = {};


for pcnt = [1, 2, 3, 4, 6, 7, 8, 9]

param_paramterization_type = params(pcnt).param_paramterization_type;    
param_polynomial_degree = params(pcnt).param_polynomial_degree;
param_num_control_points = params(pcnt).param_num_control_points;


pstr = [];
for k = 1 : length(param_polynomial_degree)
    pstr = [pstr, num2str(param_polynomial_degree(k))];
end
cstr = [];
for k = 1 : length(param_num_control_points)
    cstr = [cstr, num2str(param_num_control_points(k))];
end
if strcmp(param_paramterization_type, 'Polynomial')
    append_info = ['p', pstr];
end
if strcmp(param_paramterization_type, 'BSpline')
    append_info = ['p', pstr, '_c', cstr];
end
exp_curve_names = [exp_curve_names,  append_info];




rs_data_dir = [data_dir, '/rolling_shutter/'];

old_dir = pwd;
cd(data_dir);
gt_poses_array = getExMatrix();
cd(old_dir);



RMSE_rot_arr = [RMSE_rot_arr;   zeros(1, length(frame_names))];
RMSE_pos_arr = [RMSE_pos_arr;   zeros(1, length(frame_names))];


for ii = 1 : length(frame_names)
    
    filename = frame_names{ii};
    
    keypoints_rollingshutter = readmatrix([rs_data_dir,  filename, '.txt'])';
    keypoints_template = readmatrix([rs_data_dir, filename, '_template.txt'])';
    
    image_rollingshutter = imread([rs_data_dir, filename, '.png']);
    
    
    RSPAPP = RollingShutterPlaneAbsolutePoseProblem;
    
    RSPAPP.param_paramterization_type = param_paramterization_type;
    RSPAPP.param_polynomial_degree = param_polynomial_degree;
    RSPAPP.param_num_control_points = param_num_control_points;
    
    [rectified_img, poses] = RSPAPP.SolveImageRectification (keypoints_rollingshutter, keypoints_template, calibration_matrix, calibration_matrix, image_rollingshutter, 'image');

    
    % rectifiy the coordinate system from left-hand to right-hand
    gt_poses = gt_poses_array(ii, :);
    for tcnt = 1 : size(gt_poses, 2)
        RT = gt_poses{tcnt};
        R = RT(1:3, 1:3);
        M=diag([1, -1, -1]);
        RT(1:3, 1:3) = M * R * M;
        gt_poses{tcnt} = RT;
    end
    
    
    [RMSE_rot, RMSE_pos] = TrajectoryError.AbsoluteTrajectoryErrorByFirst (poses, gt_poses);
    RMSE_rot_arr(end, ii) = RMSE_rot;
    RMSE_pos_arr(end, ii) = RMSE_pos;
    
    
    figure('Name', append_info, 'Position', [0, 800, 1600, 600]);
    tfig = tiledlayout(1, 2, 'TileSpacing', 'tight');
    ax1 = nexttile;
    imshow(image_rollingshutter);
    ax2 = nexttile;
    imshow(rectified_img);    
    pause(1)
    
    
end


end


close all;

fontSize2 = 8;

figure('Name', 'Image rectification', 'Position', [0, 800, 1200, 400]);
tfig = tiledlayout(1, 2, 'TileSpacing', 'compact');

ax1 = nexttile;
plot(RMSE_rot_arr', '.-', 'LineWidth', 1.5);
grid on;
ylabel(ax1, 'Rotation RMSE');
% xticks(1:12);

ax2 = nexttile;
plot(RMSE_pos_arr', '.-', 'LineWidth', 1.5);
grid on;
ylabel(ax2, 'Translation RMSE');
% xticks(1:12);

lgd = legend(exp_curve_names, 'FontSize',fontSize2,'Interpreter','latex', 'Orientation', 'Horizontal', 'box', 'off');
lgd.Layout.Tile = 'South';
exportgraphics(tfig, [paper_figure_dir, 'ScanlinePoseError_',  subdir{data_choice}, '.pdf'], 'ContentType', 'Vector');




