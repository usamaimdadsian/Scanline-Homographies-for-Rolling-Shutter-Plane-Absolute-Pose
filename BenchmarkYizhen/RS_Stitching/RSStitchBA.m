function [H_RSs, H_RSs_inverse, K] = RSStitchBA(K0, tforms, num, centerImageIdx, tracks, H_vgg)

% initilization of optimization parametes 
[options, x0, lb, ub] = optSetting(K0, tforms, num, centerImageIdx, H_vgg);

% nonlinear optimziation 
[x,resnorm,residual,exitflag,output,lambda,jacobian] = ...
lsqnonlin(@unk_f_rs_homo_cost_image_v2,x0,lb,ub,options,tracks,K0,num,centerImageIdx,tforms);

% x to H_RSs 
[R,w,K] = x0toParas(x, num, centerImageIdx,K0);

for i = 2:num
    H_GS = K * R{i} * inv(K);
    A1 = K * -R{i} * X_(w{i}) * inv(K);
    A2 = K * X_(w{1}) * R{i} * inv(K);
    scale = H_GS(3,3);
    H_RSs(i).H_GS = H_GS./scale;
    H_RSs(i).A1 = A1./scale;
    H_RSs(i).A2 = A2./scale;
    
    H_GS = K * R{i}' * inv(K);
    A1 = K * -R{i}' * X_(w{1}) * inv(K);
    A2 = K * X_(w{i}) * R{i}' * inv(K);
    scale = H_GS(3,3);
    H_RSs_inverse(i).H_GS = H_GS./scale;
    H_RSs_inverse(i).A1 = A1./scale;
    H_RSs_inverse(i).A2 = A2./scale;
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% functions we need %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [options, x0, lb, ub] = optSetting(K, tforms, num, centerImageIdx, H_vgg)
    % Start with the default options
    options = optimoptions('lsqnonlin');
    options = optimoptions(options,'Display', 'iter-detailed');
    options = optimoptions(options,'Algorithm', 'levenberg-marquardt');
    options = optimoptions(options,'MaxFunEvals',100000000);
    options = optimoptions(options,'MaxIter',50);
    options = optimoptions(options,'Tolx',1e-8);
    options = optimoptions(options,'TolFun',1e-10);
    
    % parameter setting
    v_scale = [0.3 0.3 1];  % non-linear optimization parameters
    
    % camera matrix
    K(1,1) = 900;
    K(2,2) = 900;
    f0 = K(1,1);
    cols = K(1,3)*2;
    rows = K(2,3)*2;
    
    % non-linear optimiztion parameters 
    x0 = zeros(1,6*num);   % zero for angular velocity
    
    % rotations parameters from 1:i*3-3
    for i = 1:num
        % R 
        R0 = inv(K)*tforms(i).T'*K;
        %R0 = inv(K)*H_vgg{i}*K;
        angle0 = rotm2eul(R0./det(R0));
        x0(i*3-2:i*3) = [angle0(3) angle0(2) angle0(1)];  % rotation parameters
    end
    x0(centerImageIdx*3-2:centerImageIdx*3) = [];
    % focal length 
    x0 = [x0 f0];   % focal length
    
    % low and up bonds
    R0_bond = 3.14/2; w_bond = 0.2233/rows;
    f_bond_lb = 0.5*(cols+rows)/2;
    f_bond_ub = 2*(cols+rows)/2;
    
    % rotation 
    lb = -1.*[repmat(R0_bond*ones(1,3),1,num-1)];
    ub = [repmat(R0_bond*ones(1,3),1,num-1)];
    % angualr velocity 
    lb = [lb -1.*[repmat(w_bond*ones(1,3).*v_scale,1,num)]];
    ub = [ub [repmat(w_bond*ones(1,3).*v_scale,1,num)]];
    % focal length 
    lb = [lb f_bond_lb];
    ub = [ub f_bond_ub];
end


%% x0 to caemra parameters 
function [R,w,K] = x0toParas(x0, num, centerImageIdx,K)
% rotation 
for i = 1:num-1
    R{i+1} = R_caculate(x0(i*3-2),x0(i*3-1),x0(i*3));
end
% angular velocity 
w_vector = x0(3*num-2:end-1);
for i = 1:num
    w{i} = w_vector(i*3-2:i*3);
end
w_cent =  w{centerImageIdx};
w(centerImageIdx) = [];
w = {w_cent w{:}};

f0 = x0(end);
K(1,1) = f0;
K(2,2) = f0;
end