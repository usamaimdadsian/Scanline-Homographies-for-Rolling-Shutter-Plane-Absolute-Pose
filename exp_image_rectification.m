clear all
close all
clc

% this is the destination directory to contain figures for the paper
paper_figure_dir = '../PaperDraft/figure/';


size_of_final_image = [NaN, 400];

fdir = '../CVPR_dataset/';



if ~exist(paper_figure_dir, 'dir')
       mkdir(paper_figure_dir);
end

% the dirctory for the image file
%image_file = 'data_RS_for_Fang/LiU/frame02.png';

data_options = [11, 12 13, 14, 15,    21, 22, 23, 24,     31, 32,     41, 42, 43]
data_names = {'ours', 'LiU', 'Univ'};

%data_options = [11, 12 13, 14, 15,    21, 22, 23, 24,     31, 32,     41, 42, 43]
data_options = [31, 32,     42,    41, 43]
%data_options = [41, 43];

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

calibration_template = readmatrix([calibration_name, '.txt']);
calibration_rollingshutter =  readmatrix([calibration_name, '.txt']);


% solve Image rectification problem
RSPAPP = RollingShutterPlaneAbsolutePoseProblem;
rectified_img = RSPAPP.SolveImageRectification (keypoints_rollingshutter, keypoints_template, calibration_rollingshutter, calibration_template, image_rollingshutter);




image_rollingshutter = imresize(image_rollingshutter, size_of_final_image);
image_template = imresize(image_template, size_of_final_image);
rectified_img = imresize(rectified_img, size_of_final_image);

close all;



dataName = ['image_', num2str(dtc)];

fontSize2 = 8;




figure('Name', 'Scanline Homography', 'Position', [0, 1200, 550, 200]);
tfig = tiledlayout(1, 2, 'TileSpacing','Loose');

scanlineHomographies = RSPAPP.output_scanlineHomographies;
yvalues = RSPAPP.ouput_yvalues;

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
exportgraphics(tfig, [paper_figure_dir, dataName, '_Jy_estimated_curve.pdf'], 'ContentType', 'Vector');


figure('Name', 'Key Points Matching', 'Position', [0, 500, 500, 400]);
normalized_keypoints_template = RSPAPP.output_template_pts;
template_Jxpts = RSPAPP.output_template_Jxpts;
scatter(normalized_keypoints_template(1,:), normalized_keypoints_template(2,:), 'bo'); hold on;
scatter(template_Jxpts(1,:), template_Jxpts(2,:), 'r*'); hold off;
dvnorm = norm(template_Jxpts - normalized_keypoints_template, 'fro');
title(sprintf('RMSE = %f',  sqrt(dvnorm * dvnorm/ size(template_Jxpts, 2))));
xlabel('x'); ylabel('y');    


figure('Name', 'RS Image', 'Position', [700, 500, 500, 400]);
imshow(image_rollingshutter);
imwrite(image_rollingshutter, [paper_figure_dir, dataName, '_rs_img.png']);


figure('Name', 'Template Image', 'Position', [0,  0, 500, 400]);
imshow(image_template);
imwrite(image_template, [paper_figure_dir, dataName, '_template.png']);


figure('Name', 'rectified Image',  'Position', [700, 0, 500, 400]);
imshow(rectified_img);

imwrite(rectified_img, [paper_figure_dir, dataName, '_rect_img.png']);


pause(2)

end




