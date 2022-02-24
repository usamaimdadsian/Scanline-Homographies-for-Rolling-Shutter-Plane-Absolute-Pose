classdef RollingShutterPlaneAbsolutePoseProblem < handle
    
    % the tuning pramaters of the algorithm
    properties (Access = public)

        param_polynomial_degree = [1, 1, 1, 2, 2]
        
        param_paramterization_type = 'Polynomial'  %  Polynomial   or BSpline
        
        param_num_control_points = [3, 3, 3, 4, 4]
        
    end
    
    
    properties (Access = public)
        
        output_template_pts
        
        output_template_Jxpts

        ouput_yvalues
        
        output_scanlineHomographies        
        
        output_planeHomographies
        
        output_poses
        
    end
    
    
    
    methods (Access = public)
        
        
        function [rectified_img,  poses] = SolveImageRectification (this, keypoints_rollingshutter, keypoints_template, calibration_rollingshutter, calibration_template, image_rollingshutter)

            hg_rollingshutter = [keypoints_rollingshutter;  ones(1, size(keypoints_rollingshutter, 2))];
            hg_template = [keypoints_template;   ones(1, size(keypoints_template, 2))];
            
            hg_rollingshutter = inv(calibration_rollingshutter) * hg_rollingshutter;
            hg_template = inv(calibration_template) * hg_template;
            
            normalized_keypoints_rollingshutter = hg_rollingshutter([1,2], :) ./ hg_rollingshutter(3, :);
            normalized_keypoints_template = hg_template([1,2], :) ./ hg_template(3, :);

            % calculate plane-homographies in the normalized corrdinates
            JEstimate = ScanlineHomographyEstimation;
            
            JEstimate.paramterization_type = this.param_paramterization_type;
            JEstimate.polynomial_degree = this.param_polynomial_degree;            
            JEstimate.num_control_points = this.param_num_control_points;
            
            % it seems that using SVD does not improve anything
            JEstimate.use_svd_for_numerical_accuracy = false;
            
            JEstimate.EstimateScanlineHomography(normalized_keypoints_rollingshutter, normalized_keypoints_template);
            
            this.output_template_pts = normalized_keypoints_template;            
            this.output_template_Jxpts = JEstimate.WarpRS2Template(normalized_keypoints_rollingshutter);

            if exist('image_rollingshutter', 'var')
                
                [rs_y_size, rs_x_size, rs_channels] = size(image_rollingshutter);
                
                % normalize y-coordinates for each scanline
                yvalues = 1 : rs_y_size;
                hg_all_tests = [zeros(1, rs_y_size);   yvalues;  ones(1, rs_y_size)];
                hg_all_tests = inv(calibration_rollingshutter)  * hg_all_tests;
                yvalues = hg_all_tests(2, :) ./ hg_all_tests(3, :);
                
                scanlineHomographies = JEstimate.GetScanlineHomography(yvalues);
                [poses, plane_homographies] = FundamentalHomographyEquation.GetScanlinePoses (yvalues, scanlineHomographies, 100,  'image');
            
                % recalculate the un-normalized plane-homography from Euclidean plane-homographies
                for ii = 1 : length(plane_homographies)
                    EuclideanH = plane_homographies{ii};
                    plane_homographies{ii} = calibration_rollingshutter * EuclideanH * inv(calibration_template);
                end

                % use the median points of all correspondences as anchor
                anchor_idx = round(median(keypoints_rollingshutter(2, :)));
                rectified_img = ImageRectification.rectifyRSImage (image_rollingshutter, plane_homographies, anchor_idx);
                
                this.ouput_yvalues = yvalues;
                
                this.output_scanlineHomographies = scanlineHomographies;
                
                this.output_planeHomographies = plane_homographies;
                
                this.output_poses = poses;                
                
            end
            
        end
        
        
        
    end
end

