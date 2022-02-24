classdef ScanlineHomograph
    
    methods (Static, Access = public)
        
        function rota_trans = GetScanlinePose (y, scanlineHomograph)

            % the scanline pose is G_dash
            
             % M = G_dash * W
            M =  ScanlineHomograph.GetM (y, scanlineHomograph);
            
            % Recover G_dash from G_dash * W
            [R1_y, R2_y, T_y] = ScanlineHomograph.DecomposeM (M);
            
            % The last column of the rotation is determined by the first two columns
            R3_y = [-R2_y(2)*R1_y(3);
                           R2_y(1)*R1_y(3);
                           R1_y(1)/R2_y(2);];
                       
            rota_trans = [R1_y, R2_y, R3_y, T_y];
            
            % below is for debugging
            R = rota_trans(1:3, 1:3);
            if (det(R) < 0 || norm(R*R' - eye(3), 'fro') > 1e-14)
                fprintf(2, 'Pose decomposition error @ScanlineHomograph');
                fprintf(2, '\tnorm(R*R^t - eye(3)) = %.20f && det(R) = %f\n', norm(R*R' - eye(3), 'fro'), det(R));
            end
            
        end
        
        
        function rota_trans = DisambiguatePose (rota_trans, Point3D)
            point_cam = rota_trans * [Point3D(1); 0; 0; 1];
            if (point_cam(3) < 0)
                rota_trans(:, [1, 3, 4]) = - rota_trans(:, [1, 3, 4]);
                % below is for debugging
                R = rota_trans(1:3, 1:3);
                if (det(R) < 0 || norm(R*R' - eye(3), 'fro') > 1e-14)
                    fprintf(2, 'Pose decomposition error @ScanlineHomograph');
                    fprintf(2, '\tnorm(R*R^t - eye(3)) = %.20f && det(R) = %f\n', norm(R*R' - eye(3), 'fro'), det(R));
                end
            end
        end
        
        
    end

    


    
    
    methods (Static, Access = private)
        
        
        % M = G_dash * W, defined up to scale
        % G_dash is the scaline pose
        function M =  GetM (y, scanlineHomograph)

            % the common selection matrix
            W = [1, 0; 0, 0; 0, 1];
            
            % this transformation is applied to image points, 
            % this transformation vertically lifts the scanline y in image to the x-axis in image
            inv_Ddash_y = [1, 0, 0; 0, 1, y; 0, 0, 1];
            
            % model scanline
            % this is the line in the planar model (plane in 3D space, where we set z=0), corresponding to the scanlineHomograph(y)
            modelLine = cross(scanlineHomograph(:,1), scanlineHomograph(:,2)); 
            
            % transformation that brings the model-scanline (with aribitrary slope and intercept) to the (model space) line Y = 0 (horizontal line at y = 0, i.e., the x-axis)
            rota_trans = ScanlineHomograph.RelativeRotationToXAxis(modelLine);
            
            % 1D homograph, along Y = 0.
            % this maps a retina point [x, y] on scaline y, to a point [x', 0] on Y=0 in the model coordinate frame
            h22 = W' * rota_trans * scanlineHomograph;
            
            M = inv_Ddash_y * W * pinv(h22);
            
        end
        
        
        % R1_y and T_y are defined up to a common sign
        % R2_y are defined up to another sign
        % that is:  R1_y = w * R1_y;
        %               R2_y = s * R2_y;
        %                  T_y = w * T_y;
        % where,  w = +1 or -1;   s = +1 or -1;
        function [R1_y, R2_y, T_y] = DecomposeM (M)
             
            nM1 = norm(M(:, 1), 'fro');
            
            M = M/nM1;
            
            R1_y = M(:, 1);
            
            T_y = M(:, 2);
            
            R2_y = [-M(2,1); M(1,1); 0] ./ sqrt(1-M(3,1)*M(3,1));
            
        end
        
        
        % for a line L with points p, transforming the points of L by a transformation H,
        % that is: p ---> p' = H * p
        % then the line L is transformed to L' = inv(H)^T * L
        % it can be shown that <L', p'> = <L, p> = 0, that means, if p is on L, then p' is on L'
        % Here we require H to be a rigid transformation.
        % Let  invH = [R, t; 0, 0, 1] = [cos,  -sin,  t1; 
        %                                                sin,  cos,  t2;
        %                                                   0,      0,   1; ]
        % so  invH^T = [ cos,  sin,  0;
        %                         -sin,  cos,  0; 
        %                           t1,    t2,  1; ]
        % Let the source line be.  L =[x; y; z].  Then invH^T * L is
        %   x * cos  +  y * sin;
        %  -x * sin  +  y * cos;
        %   x * t1  +  y * t2  +  z;
        % The target line is  L' = [0, c, 0]. This gives us two equations
        %  x * cos  +  y * sin = 0      --------------   tan = -x/y
        %  x * t1  +  y * t2  +  z = 0  ------------- point [t1, t2, 1] is one the line L
        % if y != 0,  we let t1 = 0,  t2 = -z/y
        % if x != 0,  we let t2 = 0,  t1 = -z/x
        function T = RelativeRotationToXAxis(modelLine)
            if (abs(modelLine(2)) > 1e-14)
                % the case y ! = 0
                theta = atan(- modelLine(1)/modelLine(2));
                t1 = 0;  t2 = -modelLine(3)/modelLine(2);
            else
                % the case y = 0. the constraints become:
                % x * cos = 0
                % x * t1 + z = 0
                theta = pi/2;
                t2 = 0;  t1 = -modelLine(3)/modelLine(1);
            end
            c = cos(theta);  s = sin(theta);
            R = [c, -s;  s, c];  T = [t1; t2];
%             invT = [R, T; 0, 0, 1];
            T = [R',  -R'*T; 0, 0, 1];
        end
        

        
    end
end

