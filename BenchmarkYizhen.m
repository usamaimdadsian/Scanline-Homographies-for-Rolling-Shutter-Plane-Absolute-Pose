classdef BenchmarkYizhen < handle
 
    
    methods (Access = public)
  

        %%%%%%  Wrapper for Yizhen's method:: %%%%
        % PARAMETERS:
        % keypoints_rollingshutter: [2 x N]
        % keypoints_template: [2 x N]
        % calibration_rollingshutter: [3 x 3]
        % ---->calibration_template: [3 x 3]   - ASSUMPTION1:
        % ----> calibration_rollingshutter = calibration_template, this parameter has been DEPRECATED
        % image_rollingshutter: RS image - ASSUMPTION2: NO radial distortion in image
        % template_type: image (object type template NOT applicable)
        
        % OUTPUT:
        % rectified_img - images        
        % poses -  for each scanline, format [N x 3]
        % landmarkHomography - [N x 1] cell array of landmark homographies
        % status: boolean indicating is rectification completed or not


        function [rectified_img,  poses] = SolveImageRectification (this, keypoints_rollingshutter, keypoints_template, calibration_rollingshutter, calibration_template, image_rollingshutter, template_type)
            if ~exist('template_type', 'var') || (exist('template_type', 'var') && isempty(template_type))
                template_type = 'image';    %% image = global-shutter-image,     object = Euclidean object
            end
            
            addpath(genpath('./BenchmarkYizhen'));
            
            close all;

            poses = cell(0);
            
            rectified_img = TwoViewRectificationProcedure (keypoints_rollingshutter, keypoints_template, image_rollingshutter, [], calibration_rollingshutter);
            
            [rsy, rsx, ch] = size(image_rollingshutter);

            rectified_img = imresize(rectified_img, [rsy, rsx]);

            close all;

        end




        function [rectified_landmarks,  poses] = SolveLandmarkRectification (this, keypoints_rollingshutter, keypoints_template, calibration_rollingshutter, calibration_template, landmarks_rollingshutter, template_type)
            if ~exist('template_type', 'var') || (exist('template_type', 'var') && isempty(template_type))
                template_type = 'image';    %% image = global-shutter-image,     object = Euclidean object
            end

            addpath(genpath('./BenchmarkYizhen'));


        end
        
        
    end    
    
end

