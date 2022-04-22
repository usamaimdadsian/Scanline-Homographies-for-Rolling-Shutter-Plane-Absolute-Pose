clear all
close all
clc

addpath(genpath('./RollingShutterPlaneAbsolutePose'));

size_of_final_image = [NaN, 400];


fdir = '../CVPR_dataset/';
% the list of datasets
%data_options = [11, 12, 13, 14, 15,    21, 22, 23, 24,     31, 32,     41, 42, 43]
data_options = [31, 32,     41, 42, 43]


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



for pcnt = 1 : length(params)

param_paramterization_type = params(pcnt).param_paramterization_type;    
param_polynomial_degree = params(pcnt).param_polynomial_degree;
param_num_control_points = params(pcnt).param_num_control_points;


for dtc = data_options

switch dtc

    case 11
        
        template_name = [fdir, 'Ours/data1/template'];
        rollingshutter_name = [fdir, 'Ours/data1/frame258_m'];
        
    case 12
        
        template_name = [fdir, 'Ours/data2/template'];
        rollingshutter_name = [fdir, 'Ours/data2/frame285_m'];        
        
    case 13

        template_name = [fdir, 'Ours/data3/template'];
        rollingshutter_name = [fdir, 'Ours/data3/frame534_m'];                
        
    case 14
 
         template_name = [fdir, 'Ours/data4/template'];
        rollingshutter_name = [fdir, 'Ours/data4/frame562_m'];               
        
    case 15
        
        template_name = [fdir, 'Ours/data5/template'];
        rollingshutter_name = [fdir, 'Ours/data5/frame1037_m'];                
        
    case 21
        
        template_name = [fdir, 'Simulation/data1/image1_template'];
        rollingshutter_name = [fdir, 'Simulation/data1/RSImg'];        
        
    case 22
        
        template_name = [fdir, 'Simulation/data2/image1_template'];
        rollingshutter_name = [fdir, 'Simulation/data2/RSImg'];               
        
    case 23

        template_name = [fdir, 'Simulation/data3/image2_template'];
        rollingshutter_name = [fdir, 'Simulation/data3/RSImg'];               
        
    case 24
        
        template_name = [fdir, 'Simulation/data4/image2_template'];
        rollingshutter_name = [fdir, 'Simulation/data4/RSImg'];               
        
    case 31
        
        template_name = [fdir, 'UnivAdelaide/data1/template'];
        rollingshutter_name = [fdir, 'UnivAdelaide/data1/RS_0036'];               
        calibration_name = [fdir, 'UnivAdelaide/CalibrationMatrix'];
            
    case 32
        
        template_name = [fdir, 'UnivAdelaide/data2/template'];
        rollingshutter_name = [fdir, 'UnivAdelaide/data2/RS_0043'];               
        calibration_name = [fdir, 'UnivAdelaide/CalibrationMatrix'];
        
    case 41        
        
        template_name = [fdir, 'LiU_real/data1/frame214_template'];
        rollingshutter_name = [fdir, 'LiU_real/data1/frame249'];               
        calibration_name = [fdir, 'LiU_real/CalibrationMatrix_iPhone4'];
%           calibration_name = [fdir, 'LiU_real/CalibrationMatrix_canonS95'];
        
    case 42
        
        template_name = [fdir, 'LiU_real/data2/frame693_template'];
        rollingshutter_name = [fdir, 'LiU_real/data2/frame700'];                       
        calibration_name = [fdir, 'LiU_real/CalibrationMatrix_iPhone4'];
        
    case 43

        template_name = [fdir, 'LiU_real/data3/frame108_template'];
        rollingshutter_name = [fdir, 'LiU_real/data3/frame418'];               
        calibration_name = [fdir, 'LiU_real/CalibrationMatrix_iPhone4'];
        
end

keypoints_rollingshutter = readmatrix([rollingshutter_name, '.txt'])';
keypoints_template = readmatrix([template_name, '.txt'])';

image_rollingshutter = imread([rollingshutter_name, '.png']);
image_template = imread([template_name, '.png']);

calibration_matrix = readmatrix([calibration_name, '.txt']);

% solve Image rectification problem
RSPAPP = RollingShutterPlaneAbsolutePoseProblem;

RSPAPP.param_paramterization_type = param_paramterization_type;
RSPAPP.param_polynomial_degree = param_polynomial_degree;
RSPAPP.param_num_control_points = param_num_control_points;

rectified_img = RSPAPP.SolveImageRectification (keypoints_rollingshutter, keypoints_template, calibration_matrix, calibration_matrix, image_rollingshutter);


if (pcnt == 1 && 1)
    % Benchmark: Yizhen's method
    RSPAPP_Yizhen = BenchmarkYizhen;
    rectified_img_yizhen = RSPAPP_Yizhen.SolveImageRectification (keypoints_rollingshutter, keypoints_template, calibration_matrix, calibration_matrix, image_rollingshutter, 'image');
end



close all;

fontSize2 = 8;

dataName = ['image_', num2str(dtc)];



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





normalized_keypoints_template = RSPAPP.output_KeyPoints.templatePointsNormalized;
template_Jxpts = RSPAPP.output_KeyPoints.Jx_rollingshutterPoints;

figure('Name', 'Key Points Matching', 'Position', [0, 500, 550, 400]);
scatter(normalized_keypoints_template(1,:), normalized_keypoints_template(2,:), 'bo'); hold on;
scatter(template_Jxpts(1,:), template_Jxpts(2,:), 'r*'); hold off;
dvnorm = norm(template_Jxpts - normalized_keypoints_template, 'fro');
title(sprintf('RMSE = %f',  sqrt(dvnorm * dvnorm/ size(template_Jxpts, 2))));
xlabel('x'); ylabel('y');    
pause(0.1); 



scanlineHomographies = RSPAPP.output_RollingShutterImage.scanlineHomographies;
yvalues = RSPAPP.output_RollingShutterImage.yvaluesNormalized;

figure('Name', 'Scanline Homography', 'Position', [1200, 1200, 550, 300]);
tfig = tiledlayout(2, 1, 'TileSpacing','Loose');

linewith = {1.5, 1.5, 1.5, 1.5, 1.5, 1.5};
linecolor = {'#FF0000', '#00FF00', 	'#0000FF', 	'#00FFFF', '#FFFF00', 	'#000000'};

ax1 = nexttile;
Curve6 = zeros(6, length(scanlineHomographies));
for ii = 1 : length(scanlineHomographies)
    Curve6(:, ii) = reshape(scanlineHomographies{ii}, 6, 1);
end
hold on;
k1 = plot(yvalues,  Curve6(1, :), '-', 'LineWidth', linewith{1}, 'Color', linecolor{1});
k2 = plot(yvalues,  Curve6(2, :), '-', 'LineWidth', linewith{2}, 'Color', linecolor{2});
k3 = plot(yvalues,  Curve6(3, :), '-', 'LineWidth', linewith{3}, 'Color', linecolor{3});
k4 = plot(yvalues,  Curve6(4, :), '-', 'LineWidth', linewith{4}, 'Color', linecolor{4});
k5 = plot(yvalues,  Curve6(5, :), '-', 'LineWidth', linewith{5}, 'Color', linecolor{5});
k6 = plot(yvalues,  Curve6(6, :), '-', 'LineWidth', linewith{6}, 'Color', linecolor{6});
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
nk1 = plot(yvalues,  Curve6(1, :), '-', 'LineWidth', linewith{1}, 'Color', linecolor{1});
nk2 = plot(yvalues,  Curve6(2, :), '-', 'LineWidth', linewith{2}, 'Color', linecolor{2});
nk3 = plot(yvalues,  Curve6(3, :), '-', 'LineWidth', linewith{3}, 'Color', linecolor{3});
nk4 = plot(yvalues,  Curve6(4, :), '-', 'LineWidth', linewith{4}, 'Color', linecolor{4});
nk5 = plot(yvalues,  Curve6(5, :), '-', 'LineWidth', linewith{5}, 'Color', linecolor{5});
nk6 = plot(yvalues,  Curve6(6, :), '-', 'LineWidth', linewith{6}, 'Color', linecolor{6});
hold off;
box on;
xlabel('normalized $y$ values', 'Interpreter','latex');
ylabel('$\mathbf{\Gamma}(y)$', 'Interpreter','latex');
title('normalized $\mathbf{J}(y)$', 'Interpreter','latex')

lgd = legend([nk1, nk2, nk3, nk4, nk5, nk6], {'$\gamma_1$', '$\gamma_2$', '$\gamma_3$', '$\gamma_4$', '$\gamma_5$', '$\gamma_6$'}, ...
    'FontSize',fontSize2,'Interpreter','latex', 'Orientation', 'Horizontal', 'box', 'off');
lgd.Layout.Tile = 'South';
pause(0.1); 

if ( pcnt == length(params))
exportgraphics(tfig, [paper_figure_dir, dataName, '_Jy_estimated_curve_', append_info, '.pdf'], 'ContentType', 'Vector');
end


if ~isnan(size_of_final_image(1))
scale_t = size_of_final_image(1)/size(image_template, 1);
scale_r = size_of_final_image(1)/size(image_rollingshutter, 1);
end
if ~isnan(size_of_final_image(2))
scale_t = size_of_final_image(2)/size(image_template, 2);
scale_r = size_of_final_image(2)/size(image_rollingshutter, 2);
end

image_rollingshutter = imresize(image_rollingshutter, size_of_final_image);
image_template = imresize(image_template, size_of_final_image);
rectified_img = imresize(rectified_img, size_of_final_image);


close all;

rsfig = figure('Name', 'RS Image', 'Position', [700, 500, 500, size_of_final_image(2)]); ax = gca;
imshow(image_rollingshutter); pause(0.1); hold on;
scatter(ax, scale_r*keypoints_rollingshutter(1, :), scale_r*keypoints_rollingshutter(2, :), '.', 'g'); hold off; pause(0.1);
exportgraphics(rsfig, [paper_figure_dir, dataName, '_rs_img.pdf'], 'ContentType', 'vector');
% imwrite(image_rollingshutter, [paper_figure_dir, dataName, '_rs_img.png']);


ttfig = figure('Name', 'Template Image', 'Position', [0,  0, 500, size_of_final_image(2)]);  ax = gca;
imshow(image_template); pause(0.1); hold on;
scatter(ax, scale_t*keypoints_template(1, :), scale_t*keypoints_template(2, :), '.', 'g'); hold off; pause(0.1);
exportgraphics(ttfig, [paper_figure_dir, dataName, '_template.pdf'], 'ContentType', 'vector');
% imwrite(image_template, [paper_figure_dir, dataName, '_template.png']);


rectfig = figure('Name', 'rectified Image',  'Position', [700, 0, 500, size_of_final_image(2)]);
imshow(rectified_img); pause(0.1);
imwrite(rectified_img, [paper_figure_dir, dataName, '_rect_img_', append_info, '.png']);


if (pcnt == 1 && 1)
rectified_img_yizhen = imresize(rectified_img_yizhen, [size(rectified_img, 1), size(rectified_img, 2)]);    
rectfig_yizhen = figure('Name', 'rectified Image Yizhen',  'Position', [1400, 0, 500, size_of_final_image(2)]);
imshow(rectified_img_yizhen, 'Interpolation','bilinear'); pause(0.1);
imwrite(rectified_img_yizhen, [paper_figure_dir, dataName, '_rect_img_', 'Yizhen', '.png']);
end      

pause(1.0)

end


end
