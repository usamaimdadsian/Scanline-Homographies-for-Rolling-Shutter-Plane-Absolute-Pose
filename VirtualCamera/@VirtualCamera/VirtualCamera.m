classdef VirtualCamera < handle
    
    
    properties (Access = public)
        
        ImgPoints = struct('Data', [], 'Vis', [])
        
        CameraComponents = struct('IntrinsicMatrix', [], 'Rotation', [], 'Center', [])

        GT_CameraMatrix = {}
                
        GT_Points3D = []
 
    end
    


    properties (Access = private)
        
        num_cameras = 0
        
        num_points = 0
        
        cameraFOV = pi * 120/180
        
        cameraMinRange = 10
        
        cameraMaxRange = 200
        
        cameraPrinciplePoints = [0; 0]
        
        cameraFocalLength = []
        
        xPixelPerUnit = 2e5;  % 5 um = 5e-6.  PPD = 0.2e6
        
        yPixelPerUnit = 2e5;  % 5 um = 5e-6.   PPD = 0.2e6
        
    end
    
    
    methods (Access = public)
        
        
        function obj = VirtualCamera(Xrange, Yrange, camera_radius, image_noise_level)
            if (nargin == 0)
                Xrange = 1 : 40;
                Yrange = 1 : 20;
                camera_radius = 50;
                image_noise_level = 0;
            end
            obj.genData(Xrange, Yrange, camera_radius, image_noise_level)
        end
        
        
        function [DataCell, VisCell] = genData (this, Xrange, Yrange, camera_radius, image_noise_level)
           
            this.simulatePlanarSceneAndCameras (Xrange, Yrange, camera_radius, 0);

            this.simulateCameraInstrinsics();
            
            this.assembleCameraMatrix ();
            
            for ii = 1 : this.num_cameras
                
                this.ImgPoints(ii).Vis = this.getVisiblePoints (ii);
                            
                hgPts = this.GT_CameraMatrix{ii} * [this.GT_Points3D; ones(1, this.num_points)];
            
                normal_distribution =  randn(2, this.num_points);
                normal_distribution(normal_distribution>2.0) = 2.0;
                normal_distribution(normal_distribution<-2.0) = -2.0;
                normal_distribution(2, :) = 0;
                
                this.ImgPoints(ii).Data = hgPts([1, 2], :)./hgPts(3, :) + image_noise_level * normal_distribution ;          
                
                this.ImgPoints(ii).Data = this.ImgPoints(ii).Data .* this.ImgPoints(ii).Vis;
                
                DataCell{ii} = this.ImgPoints(ii).Data;
                
                VisCell{ii} = this.ImgPoints(ii).Vis;

            end
            
        end
        
    end

    
    methods(Access = public, Static = true)
        
        % [K, R, C] = decomposeCameraMatrix (P)
        % A function to decomposition camera matrix P to K * [R, -R*C]
        %    - K   Camera Instrinsics
        %    - R   Camera Rotation
        %    - C   Camera Center Position
        [K, R, C] = decomposeCameraMatrix (P)
        
        
        hFig = visualizePointCameraConfiguration (fid, CameraMatrices, Points3D)

        
    end    
    
    
    
    

    

    methods (Access = private)
       
        
        
        function simulatePlanarSceneAndCameras (this,  Xrange, Yrange, camera_radius, sigma_camera_pos)

            Zrange = 0;
            
            [X, Y, Z] = meshgrid(Xrange, Yrange, Zrange);
            
            this.GT_Points3D = [ vec(X')'; vec(Y')'; vec(Z')' ];

            angles = [30: 20: 150]*(pi/180);
            
            CY = camera_radius * cos(angles) + median(Yrange); 
            CZ = camera_radius * sin(angles) + median(Zrange);
            CX = median(Xrange) * ones(1, length(angles));
            
            CameraCenters = [CX;  CY; CZ];
            
            points_centroid = mean(this.GT_Points3D, 2);
            
            for ii = 1 : length(CameraCenters)
                
                this.CameraComponents(ii).Center = CameraCenters(:, ii);
                % simulate camera orientations by facing towards/outwards the centroid of the point-cloud
                
                principleAxis = points_centroid - this.CameraComponents(ii).Center;
                
                normal_distribution =  randn(3, 1);
                normal_distribution(normal_distribution>2.0) = 2.0;
                normal_distribution(normal_distribution<-2.0) = -2.0;
                
                principleAxis = principleAxis + sigma_camera_pos * normal_distribution;
                
%                 this.CameraComponents(ii).Rotation = VirtualCamera.cameraRotationByTwoAxes (principleAxis, [1; 0; 0]);
                
                this.CameraComponents(ii).Rotation = VirtualCamera.cameraRotationByPrincipleAxis (principleAxis);
                
            end
            
            this.num_points = length(this.GT_Points3D);
            
            this.num_cameras = length(this.CameraComponents);
            
        end
        
        
        
        function simulateCameraInstrinsics (this)
            
            % generate a sequence of camera focal lengths
            this.cameraFocalLength = 0.01 * [0.25 + 0.01*(0 : this.num_cameras-1)];
            
            px = this.cameraPrinciplePoints(1) * this.xPixelPerUnit;
            py = this.cameraPrinciplePoints(2) * this.yPixelPerUnit;
            
            for ii = 1 : this.num_cameras
            
                fx = this.cameraFocalLength(ii) * this.xPixelPerUnit;
                fy = this.cameraFocalLength(ii) * this.yPixelPerUnit;
                
                px = px + 0*randn;  py = py + 0*randn;

                K = [ fx,  0,  px;
                         0,  fy,  py;
                         0,  0,   1 ];
                
                this.CameraComponents(ii).IntrinsicMatrix = K;
                
            end
       
        end 

        
        
        
        % assemble camera matrix, combing together camera intrinsics and camera poses
        function assembleCameraMatrix (this)
            
            for ii = 1 : this.num_cameras
                
                K = this.CameraComponents(ii).IntrinsicMatrix;
                
                R = this.CameraComponents(ii).Rotation;
                
                t = this.CameraComponents(ii).Center;
                
                this.GT_CameraMatrix{ii} = K * [R, -R*t];  % K * R * [eye(3),  -t]
                
            end
            
        end
        
        
        
        function vis = getVisiblePoints (this, camera_idx)
            
            % transform points into the camera frame
            ptsCam = this.CameraComponents(camera_idx).Rotation * (this.GT_Points3D - this.CameraComponents(camera_idx).Center);
            
            % points in range [min_range, max_range]
            vis_min_range = (ptsCam(3, :) > this.cameraMinRange);
            vis_max_range = (ptsCam(3, :) < this.cameraMaxRange);
            % points in fov: angle with respect to zaxis
            vis_fov = (acos(ptsCam(3,:) ./ sqrt(sum(ptsCam .* ptsCam, 1))) < (this.cameraFOV/2)); 
            
            vis = logical (vis_min_range .* vis_max_range .* vis_fov);
            
        end

        
    end
    
    

    
    methods (Access = private, Static = true)
                
        
        function R = cameraRotationByPrincipleAxis (principleAxis)
            
            [~, ~, V ] = svd (principleAxis');
            
            R = V(:, [3,2,1])';
            
            if (R(3,:) * principleAxis < 0)
                R(3, :) = - R(3, :);
            end
            
            if (det(R) < 0)
                R(1, :) = - R(1, :);
            end
            
        end
        

        
        function R = cameraRotationByTwoAxes (principleAxis, anotherAxis)
            
            [~, ~, V ] = svd ([2*principleAxis/norm(principleAxis), anotherAxis]');
            
            R = V(:, [2,3,1])';
            
            if (R(3,:) * principleAxis < 0)
                R(3, :) = - R(3, :);
            end
            
            if (det(R) < 0)
                R(1, :) = - R(1, :);
            end
            
        end

    end
    
    
    
end







