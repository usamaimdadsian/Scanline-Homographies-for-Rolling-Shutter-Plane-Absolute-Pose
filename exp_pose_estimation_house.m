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

template_image_name = 'gt_start_frame09_template';
gs_template_image = imread([fdir, template_image_name, '.png']);

save_param_option = 2;



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


for pcnt = 1 : length(params)
    
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
    append_info = ['p', pstr, '\_c', cstr];
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
    
    close all;    
    
    filename = frame_names{ii};
    
    keypoints_rollingshutter = readmatrix([rs_data_dir,  filename, '.txt'])';
    keypoints_template = readmatrix([rs_data_dir, filename, '_template.txt'])';
    
    image_rollingshutter = imread([rs_data_dir, filename, '.png']);
    
    
    RSPAPP = RollingShutterPlaneAbsolutePoseProblem;
    
    RSPAPP.param_paramterization_type = param_paramterization_type;
    RSPAPP.param_polynomial_degree = param_polynomial_degree;
    RSPAPP.param_num_control_points = param_num_control_points;
    
    [rectified_img, poses] = RSPAPP.SolveImageRectification (keypoints_rollingshutter, keypoints_template, calibration_matrix, calibration_matrix, image_rollingshutter, 'image');

    
    if sum(pcnt == save_param_option && 1)
        % Benchmark: Yizhen's method
        RSPAPP_Yizhen = BenchmarkYizhen;
        rectified_img_yizhen = RSPAPP_Yizhen.SolveImageRectification (keypoints_rollingshutter, keypoints_template, calibration_matrix, calibration_matrix, image_rollingshutter, 'image');
    end
    
    close all;

    
    % rectifiy the coordinate system from left-hand to right-hand
    gt_poses = gt_poses_array(ii, :);
    for tcnt = 1 : size(gt_poses, 2)
        RT = gt_poses{tcnt};
        R = RT(1:3, 1:3);
        M=diag([1, -1, -1]);
        RT(1:3, 1:3) = M * R * M;
        gt_poses{tcnt} = RT;
    end
    
    
    [RMSE_rot, RMSE_pos] = internal_packages.TrajectoryError.RelativePoseError (poses, gt_poses);
    RMSE_rot_arr(end, ii) = RMSE_rot;
    RMSE_pos_arr(end, ii) = RMSE_pos;
    
    
    
    normalized_keypoints_template = RSPAPP.output_KeyPoints.templatePointsNormalized;
    template_Jxpts = RSPAPP.output_KeyPoints.Jx_rollingshutterPoints;
    
    figure('Name', 'Key Points Matching', 'Position', [0, 800, 500, 400]);
    scatter(normalized_keypoints_template(1,:), normalized_keypoints_template(2,:), 'bo'); hold on;
    scatter(template_Jxpts(1,:), template_Jxpts(2,:), 'r*'); hold off;
    dvnorm = norm(template_Jxpts - normalized_keypoints_template, 'fro');
    title(sprintf('RMSE = %f',  sqrt(dvnorm * dvnorm/ size(template_Jxpts, 2))));
    xlabel('x'); ylabel('y');
    title(filename);
    pause(0.1);
    
    
    
    figure('Name', append_info, 'Position', [0, 0, 1600, 600]);
    tfig = tiledlayout(1, 2, 'TileSpacing', 'tight');
    ax1 = nexttile;
    imshow(image_rollingshutter); pause(0.1); hold on;
    scatter(keypoints_rollingshutter(1, :), keypoints_rollingshutter(2, :), '.', 'r');
    hold off;
    ax2 = nexttile;
    imshow(rectified_img); pause(0.1);
    title(tfig, filename);
    
    
    if sum(pcnt == save_param_option)
        
        size_of_final_image = [NaN, 350];
        
        if ~isnan(size_of_final_image(1))
            scale_r = size_of_final_image(1)/size(image_rollingshutter, 1);
        end
        if ~isnan(size_of_final_image(2))
            scale_r = size_of_final_image(2)/size(image_rollingshutter, 2);
        end

        image_rollingshutter = imresize(image_rollingshutter, size_of_final_image);
        rectified_img = imresize(rectified_img, size_of_final_image);

        template_image = imresize(gs_template_image, size_of_final_image);

        dataName = ['House_', filename];




        ttfig = figure('Name', 'Template Image', 'Position', [0,  0, 500, size_of_final_image(2)]);
        imshow(template_image); pause(0.1); hold on;
        scatter(scale_r*keypoints_template(1, :), scale_r*keypoints_template(2, :), '.', 'g'); hold off;
        exportgraphics(ttfig, [paper_figure_dir, dataName, '_', template_image_name, '.pdf'], 'ContentType', 'vector');
        %imwrite(template_image, [paper_figure_dir, dataName, '_', template_image_name, '.pdf']);


        rsfig = figure('Name', 'RS Image', 'Position', [700, 800, 500, size_of_final_image(2)]);
        imshow(image_rollingshutter); pause(0.1); hold on;
        scatter(scale_r*keypoints_rollingshutter(1, :), scale_r*keypoints_rollingshutter(2, :), '.', 'g');
        hold off;
        pause(0.1);
        exportgraphics(rsfig, [paper_figure_dir, dataName, '_rs_img.pdf'], 'ContentType', 'vector');
        % imwrite(image_rollingshutter, [paper_figure_dir, dataName, '_rs_img.png']);
        
        rectfig = figure('Name', 'rectified Image',  'Position', [1200, 800, 500, size_of_final_image(2)]);
        imshow(rectified_img); pause(0.1);
        imwrite(rectified_img, [paper_figure_dir, dataName, '_rect_img_', append_info, '.png']);

    end
    
    
    if sum(pcnt == save_param_option && 1)
        rectified_img_yizhen = imresize(rectified_img_yizhen, [size(rectified_img, 1), size(rectified_img, 2)]);
        rectfig_yizhen = figure('Name', 'rectified Image Yizhen',  'Position', [1600, 800, 500, size_of_final_image(2)]);
        imshow(rectified_img_yizhen, 'Interpolation','bilinear'); pause(0.1);
        imwrite(rectified_img_yizhen, [paper_figure_dir, dataName, '_rect_img_', 'Yizhen', '.png']);
    end
    
    pause(1);
    
end


end




close all;

fontSize2 = 8;

figure('Name', 'Image rectification', 'Position', [0, 800, 1280, 300]);
tfig = tiledlayout(1, 2, 'TileSpacing', 'compact');

linespec = {'.-', '.-', '.-', '.-', '.-',     'o--', 'o--', 'o--', 'o--', 'o--'};
linewith = {1.5, 1.5, 1.5, 1.5, 1.5,     2, 2, 2, 2, 2};
linecolor = {'#FF0000', '#00FF00',  '#0000FF', '#00FFFF', '#FFFF00', ...
                      '#A2142F', '#77AC30',  '#0072BD', '#4DBEEE', '#EDB120'};


ax1 = nexttile;
hold on;
for ii = 1 : size(RMSE_rot_arr, 1)
plot(RMSE_rot_arr(ii, :), linespec{ii}, 'LineWidth', linewith{ii}, 'Color', linecolor{ii});
end
hold off; grid off; box on;
ylabel(ax1, 'Rotation RMSE');
ax1.XLim(1) = 1;
set(gca, 'YScale', 'log');


ax2 = nexttile;
hold on;
for ii = 1 : size(RMSE_pos_arr, 1)
plot(RMSE_pos_arr(ii, :), linespec{ii}, 'LineWidth', linewith{ii}, 'Color', linecolor{ii});
end
hold off; grid off; box on;
ylabel(ax2, 'Translation RMSE');
ax2.XLim(1) = 1;
set(gca, 'YScale', 'log');

lgd = legend(exp_curve_names, 'FontSize',fontSize2,'Interpreter','latex', 'Orientation', 'Horizontal', 'box', 'off');
lgd.Layout.Tile = 'South';
pause(0.1); 

exportgraphics(tfig, [paper_figure_dir, 'ScanlinePoseError_',  subdir{data_choice}, '.pdf'], 'ContentType', 'Vector');







