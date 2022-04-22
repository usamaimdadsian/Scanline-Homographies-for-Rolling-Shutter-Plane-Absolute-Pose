classdef VirtualRollingShutterCamera < handle
    
    
    properties (Access = public)
        
        ImgPoints = struct('RollingShutterAll', [], 'RollingShutterPerScanline', [], 'GlobalShutter', [])
        
        GT_RollingShutterCameraMatrix = struct('CalibrationMatrix', [],  'PerScanlinePoses', [],  'PerPointPoses', [])
        
        GT_Points3D = struct('All', [], 'PerScanline', [])
 
        param_rs_sigma_rot = 0.01;
        
        param_rs_sigma_pos = 0.01;
        
        global_shutter_anchor = 1
        
    end
    


    properties (Access = private)
                
        num_cameras = 0
        
        num_points = 0        
        
        cameraPrinciplePoints = [0; 0]
        
        cameraFocalLength = []
        
        xPixelPerUnit = 2e5;  % 5 um = 5e-6.  PPD = 0.2e6
        
        yPixelPerUnit = 2e5;  % 5 um = 5e-6.   PPD = 0.2e6
        
        CameraComponents = struct('IntrinsicMatrix', [], 'Rotation', [], 'Center', [])
        
    end
    
    
    properties (Access = private)
        
        num_of_scanlines = []  % for each image
        
    end
    
    
    methods (Access = public)
        
        
        function obj = VirtualRollingShutterCamera (Xrange, Yrange, camera_radius, image_noise_level, rs_param_rot, rs_param_pos)
            if (nargin == 0)
                Xrange = 1 : 40;
                Yrange = 1 :  5 : 20;
                camera_radius = 50;
                image_noise_level = 0;
            end
            if exist('rs_param_rot', 'var')
                obj.param_rs_sigma_rot = rs_param_rot;
            end
            if exist('rs_param_pos', 'var')
                obj.param_rs_sigma_pos = rs_param_pos;
            end
            obj.genData(Xrange, Yrange, camera_radius, image_noise_level);
        end
        
        
        function genData (this, Xrange, Yrange, camera_radius, image_noise_level)
           
            this.simulatePlanarSceneAndCameras (Xrange, Yrange, camera_radius, 0);

            this.simulateCameraInstrinsics();
            
            this.assembleRollingShutterCameraMatrix ();
                        
            for ii = 1 : this.num_cameras

                K = this.GT_RollingShutterCameraMatrix(ii).CalibrationMatrix;
                rs_poses = this.GT_RollingShutterCameraMatrix(ii).PerScanlinePoses;
                
                hgpts = K * rs_poses{this.global_shutter_anchor} * [this.GT_Points3D.All;  ones(1, size(this.GT_Points3D.All, 2));];
                this.ImgPoints(ii).GlobalShutter = hgpts([1, 2], :) ./ hgpts(3, :);
                
                ImgPts = [];
                
                PerPointPoses = cell(1, size(this.GT_Points3D.All, 2));
                tmp = 1;
                
                for jj = 1 : this.num_of_scanlines                   
      
                    Points3DPerScanline = this.GT_Points3D.PerScanline{jj};
                     
                    hgPts = K * rs_poses{jj} * [Points3DPerScanline; ones(1, size(Points3DPerScanline, 2))];
                    
                    normal_distribution =  randn(2, size(hgPts, 2));
                    normal_distribution(normal_distribution>3.0) = 3.0;
                    normal_distribution(normal_distribution<-3.0) = -3.0;
                    
                    this.ImgPoints(ii).RollingShutterPerScanline{jj} = hgPts([1, 2], :)./hgPts(3, :) + image_noise_level * normal_distribution ;
        
                    ImgPts = [ImgPts,  this.ImgPoints(ii).RollingShutterPerScanline{jj}];
                    
                    PerPointPoses(tmp: (tmp+size(Points3DPerScanline, 2)-1)) = rs_poses(jj);
                    tmp = tmp +size(Points3DPerScanline, 2);
                    
                end
                
                this.ImgPoints(ii).RollingShutterAll = ImgPts;                
                
                this.GT_RollingShutterCameraMatrix(ii).PerPointPoses = PerPointPoses;
                
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

            % the rolling-shutter direction is along Y
            this.num_of_scanlines = length(Yrange);
            
            Pts3D = [];
            for jj = 1 : this.num_of_scanlines                           
                [X, Y, Z] = meshgrid(Xrange, Yrange(jj), 0);
                this.GT_Points3D.PerScanline{jj} = [X; Y; Z];
                Pts3D = [Pts3D,   this.GT_Points3D.PerScanline{jj}];
            end            
            this.GT_Points3D.All = Pts3D;
            
            points_centroid = mean(Pts3D, 2);
            
            angles = [30: 20: 150]*(pi/180);
            
            CY = camera_radius * cos(angles) + median(Yrange); 
            CZ = camera_radius * sin(angles);
            CX = median(Xrange) * ones(1, length(angles));
            
            CameraCenters = [CX;  CY; CZ];
                                   
            for ii = 1 : length(CameraCenters)
                
                this.CameraComponents(ii).Center = CameraCenters(:, ii);
                % simulate camera orientations by facing towards/outwards the centroid of the point-cloud
                
                principleAxis = points_centroid - this.CameraComponents(ii).Center;
                
                normal_distribution =  randn(3, 1);
                normal_distribution(normal_distribution>2.0) = 2.0;
                normal_distribution(normal_distribution<-2.0) = -2.0;
                
                principleAxis = principleAxis + sigma_camera_pos * normal_distribution;
                
                this.CameraComponents(ii).Rotation = this.cameraRotationByTwoAxes (principleAxis, [1; 0; 0]);                
%                 this.CameraComponents(ii).Rotation = this.cameraRotationByPrincipleAxis (principleAxis);
                
            end
            
            this.num_points = length(this.GT_Points3D.All);
            
            this.num_cameras = length(this.CameraComponents);
            
        end
        
        
        
        function simulateCameraInstrinsics (this)
            
            % generate a sequence of camera focal lengths
            this.cameraFocalLength =0.001*ones(this.num_cameras,1);% 0.01 * [0.25 + 0.01*(0 : this.num_cameras-1)];
            
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
        function assembleRollingShutterCameraMatrix (this)
            
            for ii = 1 : this.num_cameras
                
                K = this.CameraComponents(ii).IntrinsicMatrix;                
                R = this.CameraComponents(ii).Rotation;                
                C = this.CameraComponents(ii).Center;
                
                GS_CameraPose = [R,  -R*C];  % K * R * [eye(3),  -C]
                
                MS = internal_packages.MotionSimulator(GS_CameraPose);
                
                simulation_steps = this.num_of_scanlines - 1;
                                
                MS.SimuStaticMotionGaussianNoise (this.param_rs_sigma_rot, this.param_rs_sigma_pos, simulation_steps);
                
                this.GT_RollingShutterCameraMatrix(ii).CalibrationMatrix = K;
                this.GT_RollingShutterCameraMatrix(ii).PerScanlinePoses = MS.mPoseCell;
                
                
            end
            
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







