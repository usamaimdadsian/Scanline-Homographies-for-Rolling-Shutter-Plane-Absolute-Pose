% Copyright 2022 Fang Bai <fang dot bai at yahoo dot com>
% 
% Licensed under the Apache License, Version 2.0 (the "License");
% you may not use this file except in compliance with the License.
% You may obtain a copy of the License at
% 
%     http://www.apache.org/licenses/LICENSE-2.0
% 
% Unless required by applicable law or agreed to in writing, software
% distributed under the License is distributed on an "AS IS" BASIS,
% WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the License for the specific language governing permissions and
% limitations under the License.

classdef BSpline
    
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
    
    methods (Access = public, Static = true)
        
        function NBasisVec = GetSplineBasisVector(query_t, knot_vector, polynomial_degree_p)            
            num_knots = length(knot_vector);
            % initialize the basis vector, for order 0
            NBasisVec = zeros(1, num_knots-1); idx = 0;
            for ii = 1 : num_knots-1
                if (query_t >= knot_vector(ii) && query_t < knot_vector(ii+1))
                    idx = ii; 
                    NBasisVec(idx) = 1;
                    break;
                end
            end
            if (idx == 0 || query_t == knot_vector(1) || query_t == knot_vector(end))
                fprintf(2, 'fatal error@BSpline. query_t (%f) is not inside the interval of knot_vector [%f, %f]!\n',  query_t,  knot_vector(1),  knot_vector(end)); return;
            end
            % Cox-de Boor recursion formula
            % recursive construction, for k = 1, ... order_k
            % NtmpVec ---> NkVec
            for p = 1 : polynomial_degree_p
                NtmpVec = zeros(1, length(NBasisVec)-1);
                for ii = 1 : (length(NBasisVec)-1)
                    Nc1 = NBasisVec(ii);
                    Nc2 = NBasisVec(ii+1);
                    %%% -- Begin Debug -- %%%
                    if (1)
                        if knot_vector(ii+p) == knot_vector(ii) && Nc1 ~= 0
                            fprintf(2, 'exception@BSpline!\n'); return;
                        end
                        if knot_vector(ii+p+1) == knot_vector(ii+1) && Nc2 ~= 0
                            fprintf(2, 'exception@BSpline!\n'); return;
                        end
                    end
                    %%% -- End Debug -- %%%
                    if (Nc1 > 0)
                        c1 = (query_t - knot_vector(ii)) / (knot_vector(ii+p) - knot_vector(ii));
                        Nc1 = c1 * Nc1;
                    end
                    if (Nc2 > 0)
                        c2 = (knot_vector(ii+p+1) - query_t) / (knot_vector(ii+p+1) - knot_vector(ii+1));
                        Nc2 = c2 * Nc2;
                    end
                    NtmpVec(ii) = Nc1 + Nc2;
                end
                NBasisVec = NtmpVec;                
            end            
            % verifiy the construction
            if (1)
                correctness = BSpline.CheckBasisFunctions (NBasisVec, polynomial_degree_p, idx);
                if (~correctness)
                    fprintf(2, 'fatal error@BSpline.\n'); return;
                end
            end
        end

        
        function knot_positions = create_knot_positions(param_range, polynomial_degree_p, num_control_points)
            
            %%% num_knots (m+1) = num_control_points (n+1) + polynomial_degree (p) + 1
            num_knots = num_control_points+polynomial_degree_p+1;
            
            %%% Clamped B-Spline
             %%% The first k=p+1 knots must be at the same position. So must be the last k=p+1 knots                        
            k = polynomial_degree_p+1;
            
            % interpolation range
            inter_lower_limit = min(param_range); 
            inter_upper_limit = max(param_range);

            % extrapolation range
            tuning_param = 5;
            span_range = inter_upper_limit - inter_lower_limit;
            extra_lower_limit = inter_lower_limit - tuning_param * span_range;
            extra_upper_limit = inter_upper_limit + tuning_param * span_range;
            
            num_rest_knots = num_knots - 2 * k;
            
            if (num_rest_knots > 0)
                % in this case, we use clamped B-spline
                ll = extra_lower_limit * ones(1, k);
                uu = extra_upper_limit * ones(1, k);                
                cc = linspace(inter_lower_limit, inter_upper_limit, num_rest_knots+2);
                cc = cc(2 : (end-1));
                knot_positions = [ll, cc, uu];
            else
                % in this case, we use unclamped B-spline
                knot_positions = linspace(inter_lower_limit, inter_upper_limit, num_knots);
                knot_positions(1) = extra_lower_limit;
                knot_positions(end) = extra_upper_limit;
            end
            
        end             
        
    end
    
    
    methods (Access = private, Static = true)
        
        function correctness =  CheckBasisFunctions (NBasisVec, polynomial_degree_p, query_idx)
            correctness = false;
            if abs(sum(NBasisVec) - 1) > 1e-14
                fprintf(2, 'fatal error@BSpline. @CheckBasisFunctions. << Partition of Unity. sum(NBasisVec) = 1 >> is violated!\n'); return;
            end
            if sum(NBasisVec((query_idx-polynomial_degree_p) : query_idx) > 0) ~= (polynomial_degree_p+1)
                fprintf(2, 'fatal error@BSpline. @CheckBasisFunctions. << NBasisVec[(query_idx-polynomial_degree_p) : query_idx] > 0 >> is violated!\n'); return;
            end
            correctness = true;
        end
    end
    
end

