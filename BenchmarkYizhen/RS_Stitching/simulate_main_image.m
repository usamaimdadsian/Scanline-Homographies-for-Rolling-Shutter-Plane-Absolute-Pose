%% RS camera under Ackermann motion
%% two views case
clc
clear all
close all

H = 480; % height of image in pixels
W = 640; % width of image in pixels
f = 320; % focal length in pixels
cx = W/2; % optical center
cy = H/2;
K = [ f 0 cx ; % Intrinsic camera parameter matrix
    0 f cy ;
    0 0 1 ];
t_v = 0.0000715;  % frame scanning speed per row

%% create 3d points
p3d = generate3DPlane([5 5], 5);   % plane 3D scene
p3d = R_caculate(pi/4,0,0) * p3d';
p3d = p3d';
p3d = p3d + repmat([0,0,25],size(p3d,1),1);

% plane property
n_vecotr = cross(p3d(2,:) - p3d(1,:), p3d(3,:) - p3d(11,:));
n_vecotr = n_vecotr./norm(n_vecotr);
n_vecotr = n_vecotr';
d =  -n_vecotr' * p3d(1,:)';
if d < 0
    d = -d;
    n_vecotr = -n_vecotr;
end
%% create two RS cameras under ackermann motion
R0 = R_caculate(0.3,-0.2,0.2);
t0 = [0;0;0];

w1 = [-5;-2;-7];
w2 = [-2;2;7];
d1 = [0;0;0];
d2 = [0;0;0];
% 
% w1 = [0;0;0];
% w2 = [0;0;0];


% 1st image
figure2D_1 = figure();
figure3D_1 = figure();
[p2d_1, camPose_1] = create_linear_RS_image(R_caculate(0,0,0), [0;0;0], w1, d1, p3d, figure2D_1, figure3D_1);

% 2nd image
figure2D_2 = figure();
figure3D_2 = figure();
[p2d_2, camPose_2] = create_linear_RS_image(R0, t0, w2, d2, p3d, figure2D_2, figure3D_2);

%% math test
H_GS = R0 - (t0*n_vecotr')./d;
H_GS_image = K* H_GS * inv(K);
scale = H_GS_image(3,3); %scale = 1;
H_GS_image = H_GS_image./scale;

A1 = -R0*X_(w1) + R0*d1*n_vecotr'./d   + t0*n_vecotr'*X_(w1)./d;
A2 = X_(w2)*R0 - d2*n_vecotr'./d;
A1_image = K* A1 * inv(K) ./ scale; 
A2_image = K* A2 * inv(K) ./ scale;

x1 = H_GS_image'; x2 = A1_image'; x3 = A2_image';
x = [x1(:);x2(:);x3(:)]; 
x(9) = [];

%% test estimation 
disp(['*****************************************************************************************************************'])
disp(['****************************RS Homography Estimation********************* '])

for i = 1:13
    u1 = p2d_1(i,1);
    u2 = p2d_2(i,1);
    v1 = p2d_1(i,2);
    v2 = p2d_2(i,2);
    
    a1 = [0,0,0, -u1, -v1, -1, u1*v2, v1*v2,...
        0,0,0, -u1*v1, -v1^2, -v1, v1*v2*u1, v2*v1^2, v1*v2...
        0,0,0, -u1*v2, -v1*v2, -v2, v2^2*u1, v2^2*v1, v2^2];
    
    a2 = [u1, v1, 1, 0,0,0, -u2*u1, -u2*v1,...
        u1*v1, v1^2, v1, 0,0,0, -u2*u1*v1, -u2*v1^2, -u2*v1...
        u1*v2, v1*v2, v2, 0,0,0, -u2*u1*v2, -u2*v1*v2, -u2*v2];
       
    A(i*2-1:i*2,:) = [a1;a2];
    b(i*2-1:i*2,:) = [-v2;u2];
   
end
A*x - b;

% GS image transformation 
disp('[GS-based Homography...............................................................')
group = [1 5 8 12];
%H_vgg = vgg_H_from_x_lin(p2d_1(group,:)',p2d_2(group,:)');
[H_vgg, inliers] = ransacfithomography_vgg([p2d_1(1:13,:)';ones(1,13)],[p2d_2(1:13,:)';ones(1,13)], 0.01);
for i = 1:size(p2d_1,1)
    u1 = p2d_1(i,1);
    u2 = p2d_2(i,1);
    v1 = p2d_1(i,2);
    v2 = p2d_2(i,2);
        
    estimatedP = H_vgg*[u1;v1;1];
    estimatedP = estimatedP(1:3)./estimatedP(3);
    p2d_2_vgg(i,1:2) =  estimatedP(1:2,1)';
end
figure(figure2D_2)
hold on 
plot(p2d_2_vgg(:,1),p2d_2_vgg(:,2),'or');
disp('********************************************************************************************************************')

% Non-linear optimization 
disp('[full model!] RS-based Non-linear optimization ..............................................................')
%% Start with the default options
options = optimoptions('lsqnonlin');
%% Modify options setting
options = optimoptions(options,'Display', 'iter-detailed');
options = optimoptions(options,'Algorithm', 'levenberg-marquardt');
options = optimoptions(options,'MaxFunEvals',100000000);
options = optimoptions(options,'MaxIter',100);       
options = optimoptions(options,'Tolx',1e-8);
options = optimoptions(options,'TolFun',1e-10);

x0 = zeros(1,9);   % zero for angular velocity
angle0 = rotm2eul(K*H_vgg*inv(K));
x0(1:3) = [angle0(3) angle0(2) angle0(1)];
x0(1:3) = [0.1 -0.2 0.3];
R0_bond = 3.14/2; w_bond = 10*t_v;
lb = -1.*[R0_bond*ones(1,3), w_bond*ones(1,6)];
ub = 1.*[R0_bond*ones(1,3), w_bond*ones(1,6)];
[x,resnorm,residual,exitflag,output,lambda,jacobian] = ...
lsqnonlin(@rs_homo_cost_image,x0,lb,ub,options,A,b,K);

disp('estimated H:')
H = R_caculate(x(1),x(2),x(3));
%H = H./H(3,3); 
H_image = K*H*inv(K);
scale = H_image(3,3);
H_image = H_image./scale;
disp('GT H (scaled already):')
%K*H_GS*inv(K)

w1_est = x(4:6);  A1_est = K*-H*X_(w1_est)*inv(K);  A1_est = A1_est./scale;
w2_est = x(7:9);  A2_est = K*X_(w2_est)*H*inv(K);   A2_est = A2_est./scale;

% non-linartr RS image trasformation 
for i = 1:size(p2d_1,1)
    u1 = p2d_1(i,1);
    v1 = p2d_1(i,2);
    
    %v2 = ((H_image(2,:)+A1_est(2,:)) * [u1;v1;1]) / (1 - A2_est(2,:)*[u1;v1;1]);
    v2 = getv2(u1,v1,H_image,A1_est,A2_est);
    H_RS = H_image + A1_est*v1 + A2_est*v2;
    estimatedP = H_RS*[u1;v1;1];
    estimatedP = estimatedP(1:3)./estimatedP(3);
    p2d_2_est(i,1:2) =  estimatedP(1:2,1)';
end
figure(figure2D_2)
hold on 
plot(p2d_2_est(:,1),p2d_2_est(:,2),'xm');



function v2 = getv2(u1,v1,H,A1,A2)
p = [u1;v1;1];
a = A2(3,:)*p;
b = H(3,:)*p + A1(3,:)*p*v1 - A2(2,:)*p;
c = -(H(2,:)+A1(2,:)*v1)*p;

v2_1 = (-b + sqrt(b^2 - 4*a*c))/(2*a);
v2_2 = (-b - sqrt(b^2 - 4*a*c))/(2*a);

if v2_1 < 0 | v2_1 >2000
    v2 = v2_2;
else 
    v2 = v2_1;
end

if a == 0 
    v2 = 0;
end
end