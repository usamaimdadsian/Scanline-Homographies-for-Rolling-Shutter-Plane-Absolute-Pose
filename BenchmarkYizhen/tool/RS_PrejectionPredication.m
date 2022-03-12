% -------------------------------------------------------------------------
% Author: Yizhen Lao
% Predication of RS projection: from 3D point P to 2D point m on RS image
%
% -------------------------------------------------------------------------
function p_img = RS_PrejectionPredication(P_w,K,R,T,w_x,d)

H = K(2,3)*2;                              % image height
fa = K(1,1);
fb = K(2,2);
u0 = K(1,3);
v0 = K(2,3);

%% parameters setup
%R0 = R';                                    % R_c_w inverse to R_w_c as R0
R0 = R;
%d = -R0 * d;                                % linear velocity
%T0 = (-R0 * T);                             % T_c_w inverse to T_w_c as T0
T0 = T;
R = w_x * R0;                               % median variables

%% Algebraic predication of projection 3D point
% determind projection row v;
% A*v^2 + B*v + C = 0
A = R(3,:) * P_w + d(3);
B = R0(3,:) * P_w + T0(3) - fb*R(2,:)*P_w - fb*d(2) - v0*R(3,:)*P_w - v0*d(3);
C = -fb*R0(2,:) * P_w - fb*T0(2) - v0*R0(3,:)*P_w - v0*T0(3);

% A = R(3,:) * P_w + d(3) - R(3,:)*T0;
% B = R0(3,:) * P_w - R0(3,:)*T0 - fb*R(2,:)*P_w + fb*R(2,:)*T0 -fb*d(2) - R(3,:)*P_w*v0 + R(3,:)*T0*v0 - d(3)*v0;
% C = -fb*R0(2,:) * P_w + fb*R0(2,:)*T0 - R0(3,:)*P_w*v0 + R0(3,:)*T0*v0;



% GS case or RS case 
if A == 0
    % GS case 
    v = -C/B;
else
    % RS case 
    % two possible solutions
    solution1 = (-B + sqrt(-4*A*C + B^2)) / (2*A);
    solution2 = (-B - sqrt(-4*A*C + B^2)) / (2*A);
    
    er_s1 = abs(solution1 - H/2);
    er_s2 = abs(solution2 - H/2);
    % pick the soluiton which fit image range
    if er_s1 < er_s2
        v = solution1;
    else
        v = solution2;
    end
end
%% Prepare output
 % determind projection column u;
u = fa * (R0(1,:)*P_w + R(1,:)*v*P_w + T0(1) + d(1)*v) / (R0(3,:)*P_w + R(3,:)*v*P_w + T0(3) + d(3)*v) + u0;
%u = fa * (R0(1,:)*P_w + R(1,:)*v*P_w - R0(1,:)*T0 - R(1,:)*T0*v + d(1)*v) /(R0(3,:)*P_w + R(3,:)*v*P_w - R0(3,:)*T0 - R(3,:)*T0*v + d(3)*v) + u0;
p_img = [u;v;1];