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
        
        function [rectified_img, poses, landmarkHomography, status] = SolveImageRectification(this, keypoints_rollingshutter, keypoints_template, calibration_rollingshutter, image_rollingshutter, template_image, landmarks_rollingshutter)
            addpath(genpath('./BenchmarkYizhen'));
            
            close all;
            
            yRS = YizhenRS(keypoints_template, keypoints_rollingshutter);
            
            yRS.set_template_image(template_image);
            yRS.set_RS_image(image_rollingshutter);
            yRS.set_intrinsics(calibration_rollingshutter);
            yRS.set_additional_landmarks(landmarks_rollingshutter);
            
            [status] = yRS.rectify();
            
            if(status)
                
                [rectified_img] = yRS.get_rectified_image();            
                [poses] = yRS.get_landmark_poses();
                [landmarkHomography] = yRS.get_landmark_homographies();
                
            else
                
                warning('Yizhen''s method for RS image rectification failed!');
                
            end
            
        end       
        
        
    end    
    
end

