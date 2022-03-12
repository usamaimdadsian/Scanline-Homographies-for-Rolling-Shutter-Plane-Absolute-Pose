classdef PlanarRSPnP < handle

    
    methods (Static = true)
        
        function [rectified_img,  poses] = SolveImageRectification (keypoints_rollingshutter, keypoints_template, calibration_rollingshutter, calibration_template, image_rollingshutter, template_type)
            if ~exist('template_type', 'var') || (exist('template_type', 'var') && isempty(template_type))
                template_type = 'object';    %% image = global-shutter-image,     object = Euclidean object
            end
            
            [R, t, Omega, v] = PlanarRSPnP.SolvePlanarSceneConstantVelocity (keypoints_rollingshutter, keypoints_template, calibration_rollingshutter);
            
            
        end
        
        
        
        function [rectified_landmarks,  poses] = SolveLandmarkRectification (keypoints_rollingshutter, keypoints_template, calibration_rollingshutter, calibration_template, landmarks_rollingshutter, template_type)
            if ~exist('template_type', 'var') || (exist('template_type', 'var') && isempty(template_type))
                template_type = 'object';    %% image = global-shutter-image,     object = Euclidean object
            end
            
            [R, t, Omega, v] = PlanarRSPnP.SolvePlanarSceneConstantVelocity (keypoints_rollingshutter, keypoints_template, calibration_rollingshutter);


            
        end

        

    end
    
    
    
    methods (Access = private,  Static = true)
        
        % template points [X, Y, 0], with Z = 0        
        function [R, t, Omega, v] = SolvePlanarSceneConstantVelocity (keypoints_rollingshutter, keypoints_template, calibration_rollingshutter)
            
            hg_rollingshutter = [keypoints_rollingshutter;  ones(1, size(keypoints_rollingshutter, 2))];            
            hg_rollingshutter = inv(calibration_rollingshutter) * hg_rollingshutter;            
            normalized_keypoints_rollingshutter = hg_rollingshutter([1,2], :) ./ hg_rollingshutter(3, :);
                   
            hg_template = [keypoints_template;   ones(1, size(keypoints_template, 2))];
            
            num_keypoints = size(keypoints_rollingshutter, 2);
            
            A = zeros(2*num_keypoints, 18);
            
            for ii = 1 : num_keypoints
                
                hgP = hg_template(:, ii);
                
                rsuv = normalized_keypoints_rollingshutter(:, ii);
                
                tmp = [hgP',    zeros(1,3),    -rsuv(1)*hgP',    rsuv(2)*hgP',   zeros(1,3),    -rsuv(1)*rsuv(2)*hgP';
                            zeros(1,3),     hgP',    -rsuv(2)*hgP',    zeros(1,3),   rsuv(2)*hgP',    -rsuv(2)*rsuv(2)*hgP'; ];
                
                A([2*ii-1, 2*ii], :) = tmp;        
                                        
            end
            
            [~, ~, V] = svd(A);   x = V(:, end);
            
            h1 = x(1:3);       h2 = x(4:6);      h3 = x(7:9);
            d1 = x(10:12);  d2 = x(13:15);  d3 = x(16:18);
            
            r1 = h1 ./ norm(h1);    r2 = h2 ./ norm(h2);   r3 = cross(r1, r2);   t = h3 ./ norm(h1);            
             
            w1 = d1 ./ norm(d1);   w2 = d2 ./ norm(d1);   v = d3 ./ norm(d1);
            
            R = [r1, r2, r3];
            
            omega_1 = ( w1(2)*R(1,2) - w2(2)*R(1,1) )  /  ( R(3,2)*R(1,1) - R(3,1)*R(1,2) );
            omega_2 = ( w1(1)*R(2,2) - w2(1)*R(2,1) )  /  ( R(3,1)*R(2,2) - R(3,2)*R(2,1) );
            omega_3 = ( w1(1)*R(3,2) - w2(1)*R(3,1) )  /  ( R(3,1)*R(2,2) - R(3,2)*R(2,1) );
            
            Omega = [ 0,    -omega_3,     omega_2;
                                omega_3,     0,     -omega_1;
                               -omega_2,    omega_1,     0;   ];        
        end
        
    end
    
    
end

