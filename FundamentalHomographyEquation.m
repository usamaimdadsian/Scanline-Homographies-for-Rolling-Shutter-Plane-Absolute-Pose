classdef FundamentalHomographyEquation

    methods (Static, Access = public)

        function [poses, homographies] = GetScanlinePoses (yvalues, scanlineHomographies, iters,  template_type)
            %
            % FUNCTION
            %       poses = GetScanlinePoses (yvalues, scanlineHomographies, iters,  template_type)
            % INPUT:
            %       @scanlineHomographies:  a cell containing scanline-homography matrices, i.e.,  J 3by2 matrices
            %       @yvalues:  a vector containing y-values of the scanlines
            %       @iters:  a scalar giving the maximum number of iterations for scanline pose refinement.
            %       @template_type.  string type. default = 'object'
            %                        'object':  tempalte is a Euclidean object
            %                        'image':  tempalte is a global-shutter image
            %
            % OUTPUT:
            %       @poses:  a cell containing all scanline poses
            % NOTE:
            %       @iters = 0.   all poses will be set to the the global-shutter approximation
            %
            if ~exist('iters', 'var') || (exist('iters', 'var') && isempty(iters))
                iters =  100;
            end
            if ~exist('template_type', 'var') || (exist('template_type', 'var') && isempty(template_type))
                template_flag = 1;  % the template is a Euclidean object
            else
                switch template_type
                    case 'object'                        
                        template_flag  =1;  % the template is a Euclidean object
                    case 'image'
                        template_flag = 2;  % the template is a global-shutter image
                end
            end
            nsize1 = length(yvalues);
            nsize2 = length(scanlineHomographies);
            if (nsize1 == nsize2)
                nsize = nsize1;
            else
                fprintf(2, 'Different sizes of yvalues and homographies!\n'); return;
            end
            
            H0 = FundamentalHomographyEquation.InitializeGlobalEuclideanHomography  (yvalues, scanlineHomographies);
            
            poses = cell(1, nsize);
            for ii = 1 : nsize
                y = yvalues(ii);
                J = scanlineHomographies{ii};
                J = J ./ norm(J(:,1), 'fro');
                H = FundamentalHomographyEquation.RefineEuclideanHomography (y, J, H0, iters);
                homographies{ii} = H;
                poses{ii} = FundamentalHomographyEquation.PoseFromEuclideanHomography(H, template_flag);
            end
        end
        
        

        % evalaluat the homography error as the optimal value of problem:
        % error = min_{scale} ||scale * H*J - [1, 0; 0, y; 0, 1]||
        % optimal scale is:  
        % scale = trace(H*J*N') / trace(H*J*J'*H');
        function err = Error_sHJ_N (y, J, H)
            %
            % INPUT:
            %       @y: the y-value of the scanline
            %       @J: scanline homography,  3 by 2 matrix
            %       @H: Plane homography,  3 by 3 matrix
            % OUTPUT:
            %       error = min_{scale} ||scale * H*J - [1, 0; 0, y; 0, 1]||
            K = H*J;
            N = [1, 0; 0, y; 0, 1];
            scale = trace(K*N') / trace(K*K');
            err = norm(scale*H*J - N, 'fro');
        end
        
    end
    
    
    
    
    
    
    
    methods (Static, Access = private)
        
        % flag = 1.    The template is a Euclidean object
        % flag = 2.    The template is a global-shuttle image, with its principale axis being orthogonal to the object plane
        function pose = PoseFromEuclideanHomography(H, flag)
            [R, S] = FundamentalHomographyEquation.DecomposeScaledH (H);
            switch flag
                case 1
                    pose = [R, R*S(:, 3)];
                case 2
                    pose = [R,  R*S(:, 3)-R(:,3)];
            end
        end

        
        % initialzie a global homography, i.e., global pose
        function H = InitializeGlobalEuclideanHomography  (yvalues, scanlineHomographies)
            nsize1 = length(yvalues);
            nsize2 = length(scanlineHomographies);
            if (nsize1 == nsize2)
                nsize = nsize1;
            else
                fprintf(2, 'Different sizes of yvalues and homographies!\n'); return;
            end
            
            Matr = zeros(3, 2*nsize);
            rhs = zeros(3, 2*nsize);
            
            for ii = 1 : nsize
                y = yvalues(ii);
                J = scanlineHomographies{ii};
                Matr(:, [2*ii-1, 2*ii]) = J ./ norm(J(:,1), 'fro');
                rhs(:, [2*ii-1, 2*ii]) = [1, 0;   0, y;   0, 1];
            end
            
            % scaledH = rhs * pinv(Matr)
            sH = (rhs * Matr') * pinv(Matr * Matr');
            
            [R, S] = FundamentalHomographyEquation.DecomposeScaledH (sH);
            
            H = R * S;
            
        end
        

        
        % iterative refinement to fovor the constraint
        % scale * R * S * J = N 
        % for each scanline pose
        function H = RefineEuclideanHomography (y, J, H0, iters)
            if ~exist('iters', 'var') || (exist('iters', 'var') && isempty(iters))
                iters = 20;
            end
            [sR, S] = FundamentalHomographyEquation.DecomposeScaledH (H0);
            for k = 1 : iters
                sR = FundamentalHomographyEquation.RefineRoationFromTranslation (y, J, S);
                S = FundamentalHomographyEquation.RefineTranslationFromRotation (y, J, sR);
                if norm(sR*S*J - [1, 0; 0, y; 0, 1], 'fro') < 1e-12
                    break;
                end
            end
            if (0)
                fprintf(1, '[%d]: Cost(R,S) ||sR*S*J-N||:  %f ',  k, norm(sR*S*J - [1, 0; 0, y; 0, 1], 'fro'));
                fprintf(1, ' f = min_s ||s*H*J - N||:  %f\n', FundamentalHomographyEquation.Error_sHJ_N (y, J, sR*S));
            end
            [R, S] = FundamentalHomographyEquation.DecomposeScaledH (sR * S);
            H = R * S;
        end
        
        
        % sR is the solution of the following scaled Procrustes problem
        % min_sR.  || s*R*S*J - N||_F^2     s.t.   R’*R = I3, det(R) = 1.
        function sR = RefineRoationFromTranslation (y, J, S)
            N = [1, 0;   0, y;   0, 1];
            [U, ~, V] = svd (S*J*N');
            R = V * U';
            if (det(R) < 0)
                R = V * diag([1, 1, -1]) * U';
            end
            scale = trace(R*S*J*N')/trace(S*(J*J')*S');
            sR = scale*R;
        end
        
        
        % S is the solution of the following problem
        % min_S.  || sR*S*J - N||_F^2,    where  S = eye(3) + [a; b; c-1]*[0, 0, 1]
        % Let v = [0, 0, 1] * J;  alpha = v * v';
        % S = eye(3) + (inv(sR)*N-J) * (v' * [0, 0, 1]) / alpha;
        function S = RefineTranslationFromRotation (y, J, sR)
            v = [0, 0, 1] * J;
            alpha = v * v';
            N = [1, 0;   0, y;   0, 1];
            scale = norm(sR(:,1), 'fro');
            R = sR/scale;
            invsR = R'/scale;
            S = eye(3) + (invsR*N-J) * (v' * [0, 0, 1]) / alpha;
        end
        
        
        % qr(scale * H) = R * (scale * S)
        function [R, S] = DecomposeScaledH (sH)
            
            [R, S] = qr(sH);
            
            sign1 = sign(S(1,1));
            R(:, 1) = sign1 * R(:, 1);
            S(1, :) = sign1 * S(1, :);
            
            sign2 = sign(S(2,2));
            R(:, 2) = sign2 * R(:, 2);
            S(2, :) = sign2 * S(2, :);
            
            if (det(R) < 0)
                R(:, 3) = - R(:, 3);
                S(3, :) = - S(3, :);
            end
            
            scale = (S(1,1)+S(2,2))/2;
            S = S / scale;
            S(1,1) = 1; S(2,2) = 1; S(1,2) = 0; 
             
        end
        
        
    end

end