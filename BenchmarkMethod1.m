classdef BenchmarkMethod1 < handle
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    % the tuning pramaters of the algorithm
    properties (Access = public)
        
        method_param1 = 0
        
        method_param2 = 0
        
    end
    
    
    % standard output for other apps
    properties (Access = public)
        
        output_template_pts
        
        output_template_Jxpts

        ouput_yvalues
        
        output_scanlineHomographies        
        
        output_planeHomographies
        
        output_poses
        
    end    
    
    
    methods (Access = public)
        
        % this returns the rectified image, and scanline poses
        function [rectified_img,  poses] = SolveImageRectification (this, keypoints_rollingshutter, keypoints_template, calibration_rollingshutter, calibration_template, image_rollingshutter, template_type)
            
            
        end
        
        
        
        function outputArg = method1(obj,inputArg)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            outputArg = obj.Property1 + inputArg;
        end
    end
end

