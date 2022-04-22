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
        
        cameraMaxRange = 35
        
        cameraPrinciplePoints = [0; 0]
        
        cameraFocalLength = []
        
        xPixelPerUnit = 2e5;  % 5 um = 5e-6.  PPD = 0.2e6
        
        yPixelPerUnit = 2e5;  % 5 um = 5e-6.   PPD = 0.2e6
        
    end
    
    
    methods (Access = public)
        
        function obj = VirtualCamera (n_num_points, n_num_cameras)
            if (nargin == 2)
                obj.num_points = n_num_points;
                obj.num_cameras = n_num_cameras;
                obj.genData(n_num_points, n_num_cameras, 0, 1);
            end
        end
        
        
        function [DataCell, VisCell] = genData (this, n_num_points, n_num_cameras, camera_noise_level, image_noise_level)
            
            this.simulatePointCameraConfigByCylinder (n_num_points, n_num_cameras, camera_noise_level);
            
            this.simulateCameraInstrinsics();
            
            this.assembleCameraMatrix ();
            
            for ii = 1 : this.num_cameras
                
                this.ImgPoints(ii).Vis = this.getVisiblePoints (ii);
                            
                hgPts = this.GT_CameraMatrix{ii} * [this.GT_Points3D; ones(1, this.num_points)];
            
                normal_distribution =  randn(2, this.num_points);
                normal_distribution(normal_distribution>2.0) = 2.0;
                normal_distribution(normal_distribution<-2.0) = -2.0;
                
                this.ImgPoints(ii).Data = hgPts([1, 2], :)./hgPts(3, :) + image_noise_level * normal_distribution;          
                
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
        
        function simulatePointCameraConfigByCylinder (this, n_num_points, n_num_cameras, sigma_camera_pos)
            
            ptsradius = 10; 
            camera_radius = 30;
            zscale = 10;

            % simulate point-cloud
            this.GT_Points3D = 2.0*ptsradius * rand(3, n_num_points) - ptsradius;
            
            % simulate camera centers
            [CX, CY, CZ] = cylinder(camera_radius, round(n_num_cameras)+0);
            CZ = zscale * CZ;
            
            Rim1 = [CX(1,1:end-1); CY(1,1:end-1); CZ(1,1:end-1); ];
            Rim1 = Rim1(:, 1:n_num_cameras);
            Rim2 = [CX(2,1:end-1); CY(2,1:end-1); CZ(2,1:end-1); ];
            Rim2 = Rim2(:, 1:n_num_cameras);
            
            points_centroid = mean(this.GT_Points3D, 2);
            
            for ii = 1 : n_num_cameras
                
                if (mod(ii,2))
                    this.CameraComponents(ii).Center =  Rim1(:, ii);
                else
                    this.CameraComponents(ii).Center =  Rim2(:, ii);
                end
%                 this.CameraComponents(ii).Center =  (Rim1(:, ii)+Rim2(:, ii))/2;
                
                % simulate camera orientations by facing towards/outwards the centroid of the point-cloud
                
                principleAxis = points_centroid - this.CameraComponents(ii).Center;
                
                normal_distribution =  randn(3, 1);
                normal_distribution(normal_distribution>2.0) = 2.0;
                normal_distribution(normal_distribution<-2.0) = -2.0;
                
                principleAxis = principleAxis + sigma_camera_pos * normal_distribution;
                
                this.CameraComponents(ii).Rotation = this.cameraRotationByPrincipleAxis (principleAxis);
                
            end
            
            this.num_points = n_num_points;
            
            this.num_cameras = n_num_cameras;
            
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
        
    end
    
    
    
end







