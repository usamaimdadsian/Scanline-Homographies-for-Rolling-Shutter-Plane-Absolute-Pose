classdef ScanlineHomographyEstimation < handle
    
    properties (Access = public)
        
        TemplatePoints = []
        
        RollingShutterImagePoints = []
        
        % the degree of the polynoimal for each B-spline
        polynomial_degree = [1, 1, 1, 2, 2]
        
        paramterization_type = 'Polynomial'        
        %paramterization_type = 'BSpline'
        
        % this number must be greater than polynomial_degree+1
        % thus the minimal value is:  polynomial_degree+2
        num_control_points = [3, 3, 3, 4, 4]
        
        % the flag to use SVD decomposition to improve numerical accuracy.
        use_svd_for_numerical_accuracy = false;
        
    end
    
        
    
    methods (Access = public)
        
        
        function obj = ScanlineHomographyEstimation(inputRollingShutterImagePoints, inputTemplatePoints)
            if exist('inputRollingShutterImagePoints', 'var') && exist('inputTemplatePoints', 'var')
                obj.RollingShutterImagePoints = inputRollingShutterImagePoints;
                obj.TemplatePoints = inputTemplatePoints;
            end
        end
        

        function JyCell = GetScanlineHomography (this, scanline_y_values)            
            if strcmp(this.paramterization_type, 'BSpline')
                B = this.BSplineCurve5(scanline_y_values, this.knot_positions, this.polynomial_degree);
            end
            if strcmp(this.paramterization_type, 'Polynomial')
                B = this.PloynomialCurve5(scanline_y_values, this.polynomial_degree);
            end            
            CurveArray = B * this.control_handles_estimate;
            for ii = 1 :  length(scanline_y_values)
                C5 = CurveArray( (5*ii-4):(5*ii) );
                J = [C5(1), C5(4);
                       C5(2), C5(5);
                        C5(3), 1];
                JyCell{ii} = J;
            end
        end
        
        
        function EstimateScanlineHomography (this, inputRollingShutterImagePoints, inputTemplatePoints)
            if exist('inputRollingShutterImagePoints', 'var') && exist('inputTemplatePoints', 'var')
                this.RollingShutterImagePoints = inputRollingShutterImagePoints;
                this.TemplatePoints = inputTemplatePoints;
            end            
            this.polynomial_degree = this.polynomial_degree(1 : this.num_spline_curves);
            
            if strcmp(this.paramterization_type, 'BSpline')
                this.num_control_points = this.num_control_points(1 : this.num_spline_curves);
                y_values = this.RollingShutterImagePoints(2, :);
                this.param_range = [min(y_values),  max(y_values)];
                % the strategy to assign knot positions.
                offset = 1e-8;
                buffer_param_range = [this.param_range(1)-offset, this.param_range(2)+offset];
                this.knot_positions = cell(1, this.num_spline_curves);
                for ii = 1 : this.num_spline_curves
                    this.knot_positions{ii} = BSpline.create_knot_positions(buffer_param_range, this.polynomial_degree(ii), this.num_control_points(ii));
                end
            end
            
            % estimate B-Spline control positions
            this.EstimateCurve5();
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
    
    
    
    
    
    properties (Access = private)
        
        num_spline_curves = 5
        
        knot_positions = cell(0, 0)
        
        control_handles_estimate = [];        
        
        param_range = [0, 1];
        
    end    
    
    
    methods (Access = private)
                 
        function EstimateCurve5 (this)
            if ( size(this.RollingShutterImagePoints, 2) == size(this.TemplatePoints, 2))
                num_key_points = size(this.RollingShutterImagePoints, 2);
            end
            
            if strcmp(this.paramterization_type, 'BSpline')
                dim_params = sum(this.num_control_points);
            end
            if strcmp(this.paramterization_type, 'Polynomial')
                dim_params = sum(this.polynomial_degree+1);
            end
            
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
                if strcmp(this.paramterization_type, 'BSpline')
                    basis_vec = this.BSplineCurve5 (y, this.knot_positions, this.polynomial_degree);
                end
                if strcmp(this.paramterization_type, 'Polynomial')
                    basis_vec = this.PloynomialCurve5 (y, this.polynomial_degree);
                end
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
        
    end
    
    
    
    
    
    methods (Access = private, Static = true)
        
        %%% By definition:   spline_order (k) = polynomial_degree (p) + 1
        %%%                            num_knots (m+1)
        %%%                            num_control_points (n+1)
        %%%
        %%% A basic relation:  m = n + p + 1,   or equivalently  m = n + k
        %%%
        %%% num_knots (m+1) = num_control_points (n+1) + polynomial_degree (p) + 1
        %%% num_knots (m+1) = num_control_points (n+1) + spline_order (k)
        
        %%% Clamped B-Spline
        %%% The first knot and the last knot must be of multiplicity k = p + 1
        %%% The first k=p+1 knots must be at the same position. So must be the last k=p+1 knots
        
        
        function B = BSplineCurve5 (query_vector, knot_positions, polynomial_degree_p)            
            internal_num_control_points = zeros(1, 5);
            for ii = 1 : 5
                internal_num_control_points(ii) = length(knot_positions{ii}) - 1 - polynomial_degree_p(ii);
            end
            
            num_queries = length(query_vector);
            
            B = zeros(5*num_queries, sum(internal_num_control_points));
            
            for idx = 1 : num_queries
                NBasisVec1 = BSpline.GetSplineBasisVector(query_vector(idx), knot_positions{1}, polynomial_degree_p(1));
                NBasisVec2 = BSpline.GetSplineBasisVector(query_vector(idx), knot_positions{2}, polynomial_degree_p(2));
                NBasisVec3 = BSpline.GetSplineBasisVector(query_vector(idx), knot_positions{3}, polynomial_degree_p(3));
                NBasisVec4 = BSpline.GetSplineBasisVector(query_vector(idx), knot_positions{4}, polynomial_degree_p(4));
                NBasisVec5 = BSpline.GetSplineBasisVector(query_vector(idx), knot_positions{5}, polynomial_degree_p(5));
                B((5*idx-4) : 5*idx, :) = blkdiag(NBasisVec1, NBasisVec2, NBasisVec3, NBasisVec4, NBasisVec5);
            end
            
        end
        
        
        
        function B = PloynomialCurve5 (query_vector, polynomial_degree_p)
            num_queries = length(query_vector);
            for idx = 1 : num_queries
                NBasisVec1 = ScanlineHomographyEstimation.PolynomialBasisVector(query_vector(idx), polynomial_degree_p(1));
                NBasisVec2 = ScanlineHomographyEstimation.PolynomialBasisVector(query_vector(idx), polynomial_degree_p(2));
                NBasisVec3 = ScanlineHomographyEstimation.PolynomialBasisVector(query_vector(idx), polynomial_degree_p(3));
                NBasisVec4 = ScanlineHomographyEstimation.PolynomialBasisVector(query_vector(idx), polynomial_degree_p(4));
                NBasisVec5 = ScanlineHomographyEstimation.PolynomialBasisVector(query_vector(idx), polynomial_degree_p(5));
                B((5*idx-4) : 5*idx, :) = blkdiag(NBasisVec1, NBasisVec2, NBasisVec3, NBasisVec4, NBasisVec5);
            end
        end
        
        
        function basisvec = PolynomialBasisVector(query_t, polynomial_degree_p)
            basisvec = ones(1,  polynomial_degree_p+1);
            for ii = 2 : length(basisvec)
                basisvec(ii) = basisvec(ii-1) * query_t;                
            end
        end
        
    end
    
    

    
end 