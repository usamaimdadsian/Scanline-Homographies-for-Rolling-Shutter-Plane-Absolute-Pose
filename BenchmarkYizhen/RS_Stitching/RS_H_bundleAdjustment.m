function [H_RSs,  camParaSet, K] = RS_H_bundleAdjustment(vSet, numImages, H_RSs, camParaSet, K)

disp('Bundle Adjustment Begin ..............................................................')

index = 1;
for i = 2:numImages
    for j = 1:i-1
        inliers = vSet.Connections.Matches{index};
        IP1 = vSet.Views.Points{i}(inliers(:,1), :)';
        IP2 = vSet.Views.Points{j}(inliers(:,2), :)';
        % forward directoin
        for k = 1:size(inliers,1)
            u1 = IP1(1,k);
            u2 = IP2(1,k);
            v1 = IP1(2,k);
            v2 = IP2(2,k);
            
            a1 = [0,0,0, -u1, -v1, -1, u1*v2, v1*v2,...
                0,0,0, -u1*v1, -v1^2, -v1, v1*v2*u1, v2*v1^2, v1*v2...
                0,0,0, -u1*v2, -v1*v2, -v2, v2^2*u1, v2^2*v1, v2^2];
            
            a2 = [u1, v1, 1, 0,0,0, -u2*u1, -u2*v1,...
                u1*v1, v1^2, v1, 0,0,0, -u2*u1*v1, -u2*v1^2, -u2*v1...
                u1*v2, v1*v2, v2, 0,0,0, -u2*u1*v2, -u2*v1*v2, -u2*v2];
            
            A{i,j}(k*2-1:k*2,:) = [a1;a2];
            b{i,j}(k*2-1:k*2,:) = [-v2;u2];
        end
        % inverse direction
        for k = 1:size(inliers,1)
            u2 = IP1(1,k);
            u1 = IP2(1,k);
            v2 = IP1(2,k);
            v1 = IP2(2,k);
            
            a1 = [0,0,0, -u1, -v1, -1, u1*v2, v1*v2,...
                0,0,0, -u1*v1, -v1^2, -v1, v1*v2*u1, v2*v1^2, v1*v2...
                0,0,0, -u1*v2, -v1*v2, -v2, v2^2*u1, v2^2*v1, v2^2];
            
            a2 = [u1, v1, 1, 0,0,0, -u2*u1, -u2*v1,...
                u1*v1, v1^2, v1, 0,0,0, -u2*u1*v1, -u2*v1^2, -u2*v1...
                u1*v2, v1*v2, v2, 0,0,0, -u2*u1*v2, -u2*v1*v2, -u2*v2];
            
            A{j,i}(k*2-1:k*2,:) = [a1;a2];
            b{j,i}(k*2-1:k*2,:) = [-v2;u2];
        end
        index = index + 1;
    end
end

% optimization prepare
[options, x0, lb, ub] = optSetting(K, numImages, camParaSet);

% iterative optimization.
% results in parameter x
[x,resnorm,residual,exitflag,output,lambda,jacobian] = ...
    lsqnonlin(@cost_func_BA,x0,lb,ub,options,A,b,numImages,K);

[H_RSs, camParaSet, K] = xProcess(x, numImages, K);
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% functions we need%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% optimization parameters setting
function [options, x0, lb, ub] = optSetting(K, numImages, camParaSet)
%iterative optimization paprameter vector
x0 = [];
% low bound and up bould
lb = [];
ub = [];

% parameter setting
f0 = K(1,1);
cols = K(1,3)*2;
rows = K(2,3)*2;
bound_scale = [0.3 0.3 1];  % non-linear optimization parameters
R0_bond = 3.14/2; w_bond = 0.2233/rows;
f_bond_lb = 0.5*(cols+rows)/2;
f_bond_ub = 2*(cols+rows)/2;

for i = 1:numImages
    % ith view: rotation
    x0(i*6-5:i*6-3) = camParaSet(i).R;
    % ith view: w
    x0(i*6-2:i*6) = camParaSet(i).w;
    % low bound
    lb(i*6-5:i*6) = -[R0_bond*ones(1,3), w_bond*ones(1,3).*bound_scale];
    % up bound
    ub(i*6-5:i*6) = [R0_bond*ones(1,3), w_bond*ones(1,3).*bound_scale];
end
% rotation first view should be fixed as I
x0(1:3) = [];
lb(1:3) = [];
ub(1:3) = [];
% add the focal length into x0 lb ub
x0 = [x0 f0];
lb = [lb f_bond_lb];
ub = [ub f_bond_ub];

% Start with the default options
options = optimoptions('lsqnonlin');
options = optimoptions(options,'Display', 'iter-detailed');
options = optimoptions(options,'Algorithm', 'levenberg-marquardt');
options = optimoptions(options,'MaxFunEvals',100000000);
options = optimoptions(options,'MaxIter',100);
options = optimoptions(options,'Tolx',1e-8);
options = optimoptions(options,'TolFun',1e-10);
end

%% after optimization, transfor x to camParaSet and f
%% prepareation
function [H_RSs, camParaSet,K] = xProcess(x0, numImages, K)


% prepare the 1st view
camParaSet(1).R = [0,0,0];
camParaSet(1).w = x0(1:3);
x0(1:3) = [];
H_RSs(1).H_GS = eye(3);
H_RSs(1).A1 = eye(3);
H_RSs(1).A2 = eye(3);

% focal length
f0 = x0(6*numImages-5);
K(1,1) = f0;
K(2,2) = f0;

% process the rest views
for n = 2:numImages
    camParaSet(n).R = x0(6*n-11:6*n-9);
    R = x0(6*n-11:6*n-9);
    R = R_caculate(R(1),R(2),R(3));
    
    camParaSet(n).w = x0(6*n-8:6*n-6);
    w = x0(6*n-8:6*n-6);
    w = X_(w);
    
    H_RSs(n).H_GS = K * R * inv(K);
    H_RSs(n).A1 = K * -R * X_(camParaSet(n).w) * inv(K);
    H_RSs(n).A2 = K * X_(camParaSet(1).w) * R * inv(K);
end

end