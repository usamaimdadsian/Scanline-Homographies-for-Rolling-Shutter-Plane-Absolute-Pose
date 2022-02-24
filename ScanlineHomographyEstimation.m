classdef ScanlineHomographyEstimation < handle

    properties (Access = public)

        TemplatePoints = []
        
        RollingShutterImagePoints = []
        
        % number of control handles. for clamped B-spline, set it to >= 5
        % 5,  means 1 knot in middle
        % 6,  means 2 knots in middle
        % 7,  means 3 knot in middle
        % 8,  means 4 knots in middle
        num_control_points = [2, 2, 2, 2, 2]
        
        % the degree of the polynoimal for each B-spline
        polynomial_order = [2, 2, 2, 2, 2]
 
         % produces R -> Rn mapping, previous default n = 5
        num_spline_curves = 5
        
        % the flag to use SVD decomposition to improve numerical accuracy.
        use_svd_for_numerical_accuracy = false; 
        
        
    end
    
        
    
    methods (Access = public)
        
        
        function obj = ScanlineHomographyEstimation(inputRollingShutterImagePoints, inputTemplatePoints)
            if exist('inputRollingShutterImagePoints', 'var') && exist('inputTemplatePoints', 'var')
                obj.RollingShutterImagePoints = inputRollingShutterImagePoints;
                obj.TemplatePoints = inputTemplatePoints;
                obj.EstimateScanlineHomography();
            end
        end
        

        function JyCell = GetScanlineHomography (this, scanline_y_values)
            if (this.num_spline_curves == 5)
                B = this.BSpline5(scanline_y_values, this.knot_positions, this.polynomial_order+1);
                CurveArray = B * this.control_handles_estimate;
                for ii = 1 :  length(scanline_y_values)
                    C5 = CurveArray( (5*ii-4):(5*ii) );
                    J = [C5(1), C5(4);
                          C5(2), C5(5);
                          C5(3), 1];
%                     J = J ./ norm(J(:, 1), 'fro');
                    JyCell{ii} = J;
                end
            end
            if (this.num_spline_curves == 6)
                B = this.BSpline6(scanline_y_values, this.knot_positions, this.polynomial_order+1);
                CurveArray = B * this.control_handles_estimate;
                for ii = 1 :  length(scanline_y_values)
                    C6 = CurveArray( (6*ii-5):(6*ii) );
                    J = [C6(1), C6(4);
                          C6(2), C6(5);
                          C6(3), C6(6)];
%                     J = J * sign(C6(6));
%                     J = J ./ norm(J(:, 1), 'fro');
                    JyCell{ii} = J;
                end
            end
        end
        
        
        function EstimateScanlineHomography (this, inputRollingShutterImagePoints, inputTemplatePoints)
            if this.num_spline_curves ~= 5 && this.num_spline_curves ~= 6
                fprintf(2, 'num_spline_curves must be 5 or 6 \n'); return;
            end
            if exist('inputRollingShutterImagePoints', 'var') && exist('inputTemplatePoints', 'var')
                this.RollingShutterImagePoints = inputRollingShutterImagePoints;
                this.TemplatePoints = inputTemplatePoints;
            end
            this.num_control_points = this.num_control_points(1 : this.num_spline_curves);
            this.polynomial_order = this.polynomial_order(1 : this.num_spline_curves);
            
            y_values = this.RollingShutterImagePoints(2, :);
            this.param_range = [min(y_values),  max(y_values)];
            
            % the strategy to assign knot positions.            
            offset = 1e-8;
            buffer_param_range = [this.param_range(1)-offset, this.param_range(2)+offset];
            
            this.knot_positions = cell(1, this.num_spline_curves);
            for ii = 1 : this.num_spline_curves
                this.knot_positions{ii} = BSpline.create_knot_positions(buffer_param_range, this.polynomial_order(ii), this.num_control_points(ii));
            end
            
            % estimate B-Spline control positions
            if (this.num_spline_curves == 5)
                this.EstimateBSpline5();
            end
            if (this.num_spline_curves == 6)
                this.EstimateBSpline6();
            end
        end
        
    end
    

    
    
    
    properties (Access = private)
        
        knot_positions = cell(0, 0)
        
        control_handles_estimate = [];        
        
        param_range = [0, 1];
        
    end    
    
    
    methods (Access = private)
                 
        function EstimateBSpline5 (this)
            if ( size(this.RollingShutterImagePoints, 2) == size(this.TemplatePoints, 2))
                num_key_points = size(this.RollingShutterImagePoints, 2);
            end
            
            dim_params = sum(this.num_control_points);

            if this.use_svd_for_numerical_accuracy
                JacobiMatr = zeros(2*num_key_points, dim_params);
                rhsVec = zeros(2*num_key_points, 1);
            else
                HessianMatr = zeros(dim_params, dim_params);
                RsdVec = zeros(dim_params, 1);
            end
                        
            for idx = 1 : num_key_points
                x = this.RollingShutterImagePoints(1, idx);
                y = this.RollingShutterImagePoints(2, idx);
                ux = this.TemplatePoints(1, idx);
                uy = this.TemplatePoints(2, idx);
                tmp = [x,   0,   -x*ux,   1,   0;
                            0,    x,   -x*uy,   0,   1; ];
                basis_vec = this.BSpline5(y, this.knot_positions, this.polynomial_order+1);
                A = tmp * basis_vec;
                if this.use_svd_for_numerical_accuracy
                    JacobiMatr((2*idx-1) : 2*idx, :) = A;
                    rhsVec((2*idx-1) : 2*idx) = [ux; uy];
                else
                    HessianMatr = HessianMatr + A'*A;
                    RsdVec = RsdVec + A' * [ux; uy];
                end
            end
            % the solution to linear least squares
            if this.use_svd_for_numerical_accuracy
                format long;  svd_Jacobian = svd(JacobiMatr)
                this.control_handles_estimate = pinv(JacobiMatr) * rhsVec;
            else
                format long;   eig_Hessian = eig(HessianMatr)
                this.control_handles_estimate = HessianMatr \ RsdVec;
            end            
        end
        
        
        function EstimateBSpline6 (this)
            if ( size(this.RollingShutterImagePoints, 2) == size(this.TemplatePoints, 2))
                num_key_points = size(this.RollingShutterImagePoints, 2);
            end
            
            dim_params = sum(this.num_control_points);
            
            if this.use_svd_for_numerical_accuracy
                CoefficientMatr = zeros(2*num_key_points, dim_params);
            else
                SystemMatr = zeros(dim_params, dim_params);
            end
            
            for idx = 1 : num_key_points
                x = this.RollingShutterImagePoints(1, idx);
                y = this.RollingShutterImagePoints(2, idx);
                ux = this.TemplatePoints(1, idx);
                uy = this.TemplatePoints(2, idx);
                tmp = [x,   0,   -x*ux,   1,   0,   -ux;
                            0,    x,   -x*uy,   0,   1,   -uy; ];
                basis_vec = this.BSpline6(y, this.knot_positions, this.polynomial_order+1);
                A = tmp * basis_vec;
                if this.use_svd_for_numerical_accuracy
                    CoefficientMatr ((2*idx-1) : 2*idx, :) = A;
                else
                    SystemMatr = SystemMatr + A'*A;
                end                
            end
            % the solution to Ralyleigh quotient
            if this.use_svd_for_numerical_accuracy
                [~, ~, V] = svd(CoefficientMatr);
                this.control_handles_estimate = V(:, end);
            else
                [V, D] = eig(SystemMatr + SystemMatr');
                [~, idx] = min(diag(D));
                this.control_handles_estimate = V(:, idx);
            end
        end
        
    end
    
    
    
    
    
    methods (Access = private, Static = true)
        
        %%% By definition:   spline_order (k) = polynomial_order (p) + 1
        %%%                            num_knots (m+1)
        %%%                            num_control_points (n+1)
        %%%
        %%% A basic relation:  m = n + p + 1,   or equivalently  m = n + k
        %%%
        %%% num_knots (m+1) = num_control_points (n+1) + polynomial_order (p) + 1
        %%% num_knots (m+1) = num_control_points (n+1) + spline_order (k)
        
        %%% Clamped B-Spline
        %%% The first knot and the last knot must be of multiplicity k = p + 1
        %%% The first k=p+1 knots must be at the same position. So must be the last k=p+1 knots
        
        function B = BSpline5 (query_vector, knot_positions, spline_order_k)
            if ~exist('spline_order_k', 'var')
                spline_order_k = 4 * ones(1, 5);     % k = 4 means cubic polynomails, i.e., 3+1.                
            end
            polynomial_order_p = spline_order_k - 1;
            
            internal_num_control_points = zeros(1, 5);
            for ii = 1 : 5
                internal_num_control_points(ii) = length(knot_positions{ii}) - 1 - polynomial_order_p(ii);
            end
            
            num_queries = length(query_vector);
            
            B = zeros(5*num_queries, sum(internal_num_control_points));
            
            for idx = 1 : num_queries
                NBasisVec1 = BSpline.GetSplineBasisVector(query_vector(idx), knot_positions{1}, polynomial_order_p(1));
                NBasisVec2 = BSpline.GetSplineBasisVector(query_vector(idx), knot_positions{2}, polynomial_order_p(2));
                NBasisVec3 = BSpline.GetSplineBasisVector(query_vector(idx), knot_positions{3}, polynomial_order_p(3));
                NBasisVec4 = BSpline.GetSplineBasisVector(query_vector(idx), knot_positions{4}, polynomial_order_p(4));
                NBasisVec5 = BSpline.GetSplineBasisVector(query_vector(idx), knot_positions{5}, polynomial_order_p(5));
                B((5*idx-4) : 5*idx, :) = blkdiag(NBasisVec1, NBasisVec2, NBasisVec3, NBasisVec4, NBasisVec5);
            end

        end

        
        function B = BSpline6 (query_vector, knot_positions, spline_order_k)
            if ~exist('spline_order_k', 'var')
                spline_order_k = 4 * ones(1, 6);     % k = 4 means cubic polynomails, i.e., 3+1.                
            end
            polynomial_order_p = spline_order_k - 1;
            
            internal_num_control_points = zeros(1, 6);
            for ii = 1 : 6
                internal_num_control_points(ii) = length(knot_positions{ii}) - 1 - polynomial_order_p(ii);
            end
            
            num_queries = length(query_vector);
            
            B = zeros(6*num_queries, sum(internal_num_control_points));
            
            for idx = 1 : num_queries
                NBasisVec1 = BSpline.GetSplineBasisVector(query_vector(idx), knot_positions{1}, polynomial_order_p(1));
                NBasisVec2 = BSpline.GetSplineBasisVector(query_vector(idx), knot_positions{2}, polynomial_order_p(2));
                NBasisVec3 = BSpline.GetSplineBasisVector(query_vector(idx), knot_positions{3}, polynomial_order_p(3));
                NBasisVec4 = BSpline.GetSplineBasisVector(query_vector(idx), knot_positions{4}, polynomial_order_p(4));
                NBasisVec5 = BSpline.GetSplineBasisVector(query_vector(idx), knot_positions{5}, polynomial_order_p(5));
                NBasisVec6 = BSpline.GetSplineBasisVector(query_vector(idx), knot_positions{6}, polynomial_order_p(6));
                B((6*idx-5) : 6*idx, :) = blkdiag(NBasisVec1, NBasisVec2, NBasisVec3, NBasisVec4, NBasisVec5, NBasisVec6);
            end

        end
        
    end
    
    
    methods (Access = public)
        
        function template_points_Jxpts = WarpRS2Template (this,  QueryRollingShutterPoints)
            x_values = QueryRollingShutterPoints(1, :);
            y_values = QueryRollingShutterPoints(2, :);            
            JyCell = this.GetScanlineHomography (y_values);            
            template_points_Jxpts = zeros(3, length(y_values));            
            for ii = 1 : length(y_values)
                J = JyCell{ii};
                template_points_Jxpts(:, ii) = J * [x_values(ii);  1];                
            end
            template_points_Jxpts = template_points_Jxpts([1,2], :) ./ template_points_Jxpts(3, :);
        end
    
    end
    
end 