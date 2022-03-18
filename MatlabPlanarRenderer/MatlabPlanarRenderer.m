classdef MatlabPlanarRenderer<handle
    %% Method to directly render a planar image texture using pre-defined 
    % 3D transformations with given intrinsics for a synthetic perspective
    % camera. The object 3D template has to be, necessarily, planar. The
    % provided texture gets fitted to the template [height x width]. Note
    % that the camera remains STATIONARY at origin (0,0,0) pointing towards +ve
    % Z-axis while the input transformation matrix moves the planar
    % template in R3
    %
    % USAGE: 
    %    tex = imread('licensePlates.jpg');
    %    mPR = MatlabPlanarRenderer(tex);
    %    T = [1 0 0 -150; 0 1 0 -150; 0 0 1 1100.5; 0 0 0 1];
    %    R = rotx(30)*roty(40)*rotz(25);
    %    T(1:3,1:3) = R;
    %    mPR.setTransformation(T);
    %    [img] = mPR.render();
    %    imshow(img)
    %    
    %   Point correspondences can be obtained from: mPR.getCorrespondence()
    %
    %   Also check the usage of:  setIntrinsics, toggleVerbosity,
    %   setDensity, getCorrespondence
    %
    %   Point correspondence format is described in the code comments
    
    
    properties(Access = private)
        texture_ = []; % image texture
        
        % intrinsics (set to some arbitrary, default values)
        height_ = 480; width_ = 640;
        fx_ = 400.0; fy_ = 400.0;
        cx_ = 240.0; cy_ = 320.0;
        
        pointDensity_ = 2;   % ideally \in [1,10]     
        
        T_ = eye(4,4);        
        verbose_ = false;
        
        correspondence_ = []; % template to pointcloud to GS image correspondence
        % DESCRIPTION of 'correspondence_':
        % Has three components:
        %    1) pointcloud = [N x 3] template points
        %    2) pointcloud_transformed - [N x 3] transformed template
        %         points, as viewed by the perspective camera
        %    3) rectified_image_indices - [H x W] array of indices
        %         corresponding to the 'pointcloud' and 'pointcloud_transformed'
        %         arrays above
        %
        %    E.g: If rectified_image_indices(x, y) = n;
        %            Then: a) pointcloud(n, :) \in R^3 is the object template
        %                           coordinates of that point
        %                      b) pointcloud_transformed(n,:) \in R^3 is the
        %                           coordinate of the point on the transformed 3D plane in camera coordinates
        %                      c) Obviously, (x,y) are the global-shutter image coordinates
        
    end
    
    
    methods(Access = public)
        
        function [obj] = MatlabPlanarRenderer(texture)
            %% Constructor, sets the image texture to be used
            obj.texture_ = texture;
        end
        
        function [obj] = setTransformation(obj,T)
            %% Set the [4x4] transformation matrix
            obj.T_ = T;
        end
        
        function [obj] = setDensity(obj, density)
            %% Set the density of the rendered plane, interpreted as a  multiplier on the image [height x width]
            obj.pointDensity_ = density;
        end
        
        function [obj] = toggleVerbosity(obj)
            %% Toggle between printing debug logs and not printing them
            obj.verbose_ = ~obj.verbose_;
        end
        
        function [obj] = setIntrinsics(obj, fx, fy, height, width)
            %% Set camera intrinsics
            obj.fx_ = fx;
            obj.fy_ = fy;
            obj.height_ = height;
            obj.width_= width;
            obj.cx_ = height/2;
            obj.cy_ = width/2;
        end
        
        function [img] = render(obj)
            %% Render the scene that has been already setup
            [img] = obj.renderImage();
        end
        
        function [correspondence] = getCorrespondence(obj)
            %% Returns the correspondence, detailed description given in Line: 41 - 56
            correspondence = obj.correspondence_;
        end
        
    end
    
    methods(Access = private)
        function [rendered_image] = renderImage(obj)
            %% Image rendering method
            % Given: [H,W] = object height, width. D = point density
            %
            % 1) We first create a virtual grid of size [HD, WD]
            % 2) We resize the image into size [HD, WD]
            % 3) We transform this grid in R3 by the the given
            %       transformation matrix T_
            % 4) We project every point of the grid to a 2D image,
            %      perspective projection from the method 'project'
            % 5) We build the rendered image from the transformed points
            %      and the point-to-texture correspondence obtained before
            
            gridHeight = obj.height_*obj.pointDensity_;
            gridWidth = obj.width_*obj.pointDensity_;
            
            texture = imresize(obj.texture_,[gridHeight gridWidth]);
            
            if(obj.verbose_) % debug block
                fprintf('Size of upsampled texture:\n');
                size(texture)            
                fprintf('Grid size: [%d x%d]\n',gridHeight, gridWidth);
            end
            
            pointcloud = ones(gridHeight*gridWidth,4);
            colorMap = ones(gridHeight*gridWidth,3);
            
            count = 1;
            
            for row = 1:gridHeight
                for col = 1:gridWidth
                    pointcloud(count,:) = [row col 0 1];
                    colorMap(count,:) = texture(row,col,:);
                    count = count + 1;
                end
            end
            
            pointcloud_transformed = obj.T_*pointcloud.';
            pointcloud_transformed = pointcloud_transformed.';
            
            if(obj.verbose_) % debug block
                fprintf('About to apply the transformation:\n');
                obj.T_
                pointcloud_debug = pointcloud_transformed(:,1:3);
                pcshow(pointcloud_debug); hold on;
                xlabel('X-axis'); ylabel('Y-axis'); zlabel('Z-axis');
                title('Densely sampled points on the transformed plane');
                pause(1);
                close all;
            end
            
            rendered_image=ones([obj.height_ obj.width_ 3],class(texture));
            rendered_image_indices = ones([obj.height_ obj.width_]);
            
            for ii = 1:(count-1)
                P = pointcloud_transformed(ii,:);
                [p] = obj.project(P(1), P(2), P(3));
                
                if((p(1) >= 1) && (p(2) >= 1) && (p(1) <= obj.height_) && (p(2) <= obj.width_))
                    rendered_image(round(p(1)), round(p(2)), :) = colorMap(ii,:);
                    rendered_image_indices(round(p(1)), round(p(2)), 1) = ii;
                end
                
            end
            
            correspondence{1} = pointcloud(:,1:3);
            correspondence{2} = pointcloud_transformed(:,1:3);
            correspondence{3} = rendered_image_indices;
            
            obj.correspondence_ = correspondence;
            
             if(obj.verbose_) % debug block
                imshow(rendered_image);
                pause(1);
                close all;
            end
            
            
        end
        
        function [p] = project(obj, X, Y, Z)
            %% Perspective projection, ideal model,  no radial/tangential distortion or skew
            
            x = (obj.fx_*(X/Z)) + obj.cx_;
            y = (obj.fy_*(Y/Z)) + obj.cy_;
            
            p = [x y];
            
        end
    end
end