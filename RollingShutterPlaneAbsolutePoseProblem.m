classdef RollingShutterPlaneAbsolutePoseProblem < handle
    
    % the tuning pramaters of the algorithm
    properties (Access = public)

        param_polynomial_degree = [1, 1, 1, 2, 2]
        
        param_paramterization_type = 'Polynomial'  %  Polynomial   or BSpline
        
        param_num_control_points = [3, 3, 3, 4, 4]
        
    end
    
    
    properties (Access = public)
        
        output_KeyPoints = struct('rollingshutterPointsNormalized',  [], 'templatePointsNormalized',  [], 'Jx_rollingshutterPoints', []);
        
        output_RollingShutterImage = struct('yvaluesNormalized', [], 'scanlineHomographies', [], 'planeHomographies', [],  'poses', []);
        
        output_RollingShutterLandmarks = struct('yvaluesNormalized', [], 'scanlineHomographies', [], 'planeHomographies', [],  'poses', [], ...
                                                                                    'rollingshutterLandmarksNormalized', [], 'Jx_rollingshutterLandmarks', []);
        
    end
    
    
    properties (Access = private)
        
        use_svd_for_numerical_accuracy = true
        
    end
    
    
    methods (Access = public)
        
        
        function [rectified_img,  poses] = SolveImageRectification (this, keypoints_rollingshutter, keypoints_template, calibration_rollingshutter, calibration_template, image_rollingshutter, template_type)
            if ~exist('template_type', 'var') || (exist('template_type', 'var') && isempty(template_type))
                template_type = 'image';    %% image = global-shutter-image,     object = Euclidean object
            end

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
            JEstimate.use_svd_for_numerical_accuracy = this.use_svd_for_numerical_accuracy;
            
            JEstimate.EstimateScanlineHomography(normalized_keypoints_rollingshutter, normalized_keypoints_template);
            
            this.output_KeyPoints.rollingshutterPointsNormalized = normalized_keypoints_rollingshutter;
            this.output_KeyPoints.templatePointsNormalized = normalized_keypoints_template;
            this.output_KeyPoints.Jx_rollingshutterPoints = JEstimate.WarpRS2Template(normalized_keypoints_rollingshutter);
            
            if exist('image_rollingshutter', 'var')
                
                [rs_y_size, rs_x_size, rs_channels] = size(image_rollingshutter);
                
                % normalize y-coordinates for each scanline
                yvalues = 1 : rs_y_size;
                hg_all_tests = [zeros(1, rs_y_size);   yvalues;  ones(1, rs_y_size)];
                hg_all_tests = inv(calibration_rollingshutter)  * hg_all_tests;
                yvalues = hg_all_tests(2, :) ./ hg_all_tests(3, :);
                
                scanlineHomographies = JEstimate.GetScanlineHomography(yvalues);
                [poses, plane_homographies] = FundamentalHomographyEquation.GetScanlinePoses (yvalues, scanlineHomographies, 100,  template_type);
            
                % recalculate the un-normalized plane-homography from Euclidean plane-homographies
                for ii = 1 : length(plane_homographies)
                    EuclideanH = plane_homographies{ii};
                    plane_homographies{ii} = calibration_rollingshutter * EuclideanH * inv(calibration_template);
                end

                % use the median points of all correspondences as anchor
                anchor_idx = round(median(keypoints_rollingshutter(2, :)));
                rectified_img = ImageRectification.rectifyRSImage (image_rollingshutter, plane_homographies, anchor_idx);
                
                this.output_RollingShutterImage.yvaluesNormalized = yvalues;
                this.output_RollingShutterImage.scanlineHomographies = scanlineHomographies;
                this.output_RollingShutterImage.planeHomographies = plane_homographies;
                this.output_RollingShutterImage.poses = poses;
                                                
            end
            
        end
        

        
        
        function [rectified_landmarks,  poses] = SolveLandmarkRectification (this, keypoints_rollingshutter, keypoints_template, calibration_rollingshutter, calibration_template, landmarks_rollingshutter, template_type)
            if ~exist('template_type', 'var') || (exist('template_type', 'var') && isempty(template_type))
                template_type = 'image';    %% image = global-shutter-image,     object = Euclidean object
            end
            
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
            JEstimate.use_svd_for_numerical_accuracy = this.use_svd_for_numerical_accuracy;
            
            JEstimate.EstimateScanlineHomography(normalized_keypoints_rollingshutter, normalized_keypoints_template);
            
            this.output_KeyPoints.rollingshutterPointsNormalized = normalized_keypoints_rollingshutter;
            this.output_KeyPoints.templatePointsNormalized = normalized_keypoints_template;
            this.output_KeyPoints.Jx_rollingshutterPoints = JEstimate.WarpRS2Template(normalized_keypoints_rollingshutter);            
            
            if exist('landmarks_rollingshutter', 'var')
                
                rs_pts_size = size(landmarks_rollingshutter, 2);
                
                % normalize y-coordinates for each scanline
                hg_all_tests = [landmarks_rollingshutter;  ones(1, rs_pts_size)];
                hg_all_tests = inv(calibration_rollingshutter)  * hg_all_tests;
                normalized_landmarks_rollingshutter = hg_all_tests([1, 2], :) ./ hg_all_tests(3, :);
                yvalues = normalized_landmarks_rollingshutter(2, :);
                
                scanlineHomographies = JEstimate.GetScanlineHomography(yvalues);
                [poses, plane_homographies] = FundamentalHomographyEquation.GetScanlinePoses (yvalues, scanlineHomographies, 100,  template_type);
                
                % recalculate the un-normalized plane-homography from Euclidean plane-homographies
                for ii = 1 : length(plane_homographies)
                    EuclideanH = plane_homographies{ii};
                    plane_homographies{ii} = calibration_rollingshutter * EuclideanH * inv(calibration_template);
                end
                
                anchor_idx = 1;
                rectified_landmarks = zeros(2, rs_pts_size);
                H0 = plane_homographies{anchor_idx};
                for ii = 1 : rs_pts_size
                    pt = landmarks_rollingshutter(:, ii);
                    H = plane_homographies{ii};
                    hgq = H0 * inv(H) * [pt; 1];
                    rectified_landmarks(:, ii) = hgq([1,2]) ./ hgq(3);
                end
                
                this.output_RollingShutterLandmarks.yvaluesNormalized = yvalues;
                this.output_RollingShutterLandmarks.scanlineHomographies = scanlineHomographies;
                this.output_RollingShutterLandmarks.planeHomographies = plane_homographies;
                this.output_RollingShutterLandmarks.poses = poses;
                                
                this.output_RollingShutterLandmarks.rollingshutterLandmarksNormalized = normalized_landmarks_rollingshutter;
                this.output_RollingShutterLandmarks.Jx_rollingshutterLandmarks = JEstimate.WarpRS2Template(normalized_landmarks_rollingshutter);            
                
            end
            
        end        
        
        
        
    end
end

