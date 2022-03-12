% Real RS images stitching  
% Author: Yizhen LAO
% Two views example 
clc 
close all
clear all



path = '.\data\syntheic\';
image1 = imread([path,'frame10.png']);
image2 = imread([path,'frame05.png']);

path = '.\data\limos\';
image1 = imread([path,'3.jpg']);
image2 = imread([path,'4.jpg']);

figure(), imshow(image1); title('Original input 1')
figure(), imshow(image2); title('Original input 2')


rows=size(image1,1); 
cols=size(image1,2);
f = 0.9 * (rows+cols)/2;
K = [f 0 cols/2;
    0 f rows/2;
    0 0 1];

%% features detection 
I1 = rgb2gray(image1);
I2 = rgb2gray(image2);
points1 = detectSURFFeatures(I1);
points2 = detectSURFFeatures(I2);
[f1,vpts1] = extractFeatures(I1,points1);
[f2,vpts2] = extractFeatures(I2,points2);
indexPairs = matchFeatures(f1,f2) ;
matchedPoints1 = vpts1(indexPairs(:,1));
matchedPoints2 = vpts2(indexPairs(:,2));
IP1 = double(matchedPoints1.Location');
IP2 = double(matchedPoints2.Location');

%% GS image transformation 
tic
disp('[GS-based Homography...............................................................')
[H_vgg, inliers] = ransacfithomography_vgg([IP1;ones(1,matchedPoints1.Count)],[IP2;ones(1,matchedPoints2.Count)], 0.01);
[image1_warp, R1] = imwarp(image1,projective2d(H_vgg'));
R2 = imref2d(size(image2));
%image_merge = imfuse(image1_warp, R1, image2, R2, 'method', 'blend');
mode = 'Average';

% move B in respect to the top-left corner of A
offsetW = round(-R1.XWorldLimits(1));
offsetH = round(-R1.YWorldLimits(1));
image_merge = blendMode(image1_warp, image2, mode, offsetW, offsetH);

figure(), imshow(image_merge); title('GS Stitching')
toc
disp('********************************************************************************************************************')


%% RS image transformation 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic
%% inverse H A1 A2 of RS 
% matching points matrix prepare 
for i = 1:size(inliers,2)
    u1 = IP2(1,inliers(i));
    u2 = IP1(1,inliers(i));
    v1 = IP2(2,inliers(i));
    v2 = IP1(2,inliers(i));
    
    a1 = [0,0,0, -u1, -v1, -1, u1*v2, v1*v2,...
        0,0,0, -u1*v1, -v1^2, -v1, v1*v2*u1, v2*v1^2, v1*v2...
        0,0,0, -u1*v2, -v1*v2, -v2, v2^2*u1, v2^2*v1, v2^2];
    
    a2 = [u1, v1, 1, 0,0,0, -u2*u1, -u2*v1,...
        u1*v1, v1^2, v1, 0,0,0, -u2*u1*v1, -u2*v1^2, -u2*v1...
        u1*v2, v1*v2, v2, 0,0,0, -u2*u1*v2, -u2*v1*v2, -u2*v2];
       
    A(i*2-1:i*2,:) = [a1;a2];
    b(i*2-1:i*2,:) = [-v2;u2];
end
% Non-linear optimization 
disp('[full model!] RS-based Non-linear optimization ..............................................................')
% Start with the default options
options = optimoptions('lsqnonlin');
% Modify options setting
options = optimoptions(options,'Display', 'iter-detailed');
options = optimoptions(options,'Algorithm', 'levenberg-marquardt');
options = optimoptions(options,'MaxFunEvals',100000000);
options = optimoptions(options,'MaxIter',100);       
options = optimoptions(options,'Tolx',1e-8);
options = optimoptions(options,'TolFun',1e-10);

% non-linear optimiztion start 
x0 = zeros(1,10);   % zero for angular velocity
% R0_vgg = inv(K)*H_vgg*K;
% angle0 = rotm2eul(R0_vgg./det(R0_vgg));
% x0(1:3) = [angle0(3) angle0(2) angle0(1)];
%x0(1:3) = [0.1 -0.2 0.3];
R0_bond = 3.14/2; w_bond = 0.2233/size(image1,1); v_scale = [0.3 0.3 1 0.3 0.3 1]; 
f_bond_lb = 0.9*(cols+rows)/2;
f_bond_ub = 2*(cols+rows)/2;
lb = -1.*[R0_bond*ones(1,3), w_bond*ones(1,6).*v_scale]; lb = [lb f_bond_lb];
ub = 1.*[R0_bond*ones(1,3), w_bond*ones(1,6).*v_scale];  ub = [ub f_bond_ub];
[x,resnorm,residual,exitflag,output,lambda,jacobian] = ...
lsqnonlin(@unk_f_rs_homo_cost_image,x0,lb,ub,options,A,b,K);

% H A1 A2 caculation 
disp('estimated H:')
H = R_caculate(x(1),x(2),x(3));
f0 = x(10);
K(1,1)=f0;
K(2,2)=f0;
H_image = K*H*inv(K);
scale = H_image(3,3);
H_image_inverse = H_image./scale;

w1_est_inverse = x(4:6);  A1_est = K*-H*X_(w1_est_inverse)*inv(K);  A1_est_inverse = A1_est./scale;
w2_est_inverse = x(7:9);  A2_est = K*X_(w2_est_inverse)*H*inv(K);   A2_est_inverse = A2_est./scale;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i = 1:size(inliers,2)
    u1 = IP1(1,inliers(i));
    u2 = IP2(1,inliers(i));
    v1 = IP1(2,inliers(i));
    v2 = IP2(2,inliers(i));
    
    a1 = [0,0,0, -u1, -v1, -1, u1*v2, v1*v2,...
        0,0,0, -u1*v1, -v1^2, -v1, v1*v2*u1, v2*v1^2, v1*v2...
        0,0,0, -u1*v2, -v1*v2, -v2, v2^2*u1, v2^2*v1, v2^2];
    
    a2 = [u1, v1, 1, 0,0,0, -u2*u1, -u2*v1,...
        u1*v1, v1^2, v1, 0,0,0, -u2*u1*v1, -u2*v1^2, -u2*v1...
        u1*v2, v1*v2, v2, 0,0,0, -u2*u1*v2, -u2*v1*v2, -u2*v2];
       
    A(i*2-1:i*2,:) = [a1;a2];
    b(i*2-1:i*2,:) = [-v2;u2];
end

% Non-linear optimization 
disp('[full model!] RS-based Non-linear optimization ..............................................................')
% Start with the default options
options = optimoptions('lsqnonlin');
% Modify options setting
options = optimoptions(options,'Display', 'iter-detailed');
options = optimoptions(options,'Algorithm', 'levenberg-marquardt');
options = optimoptions(options,'MaxFunEvals',100000000);
options = optimoptions(options,'MaxIter',100);       
options = optimoptions(options,'Tolx',1e-8);
options = optimoptions(options,'TolFun',1e-10);

% non-linear optimiztion start 
x0 = zeros(1,10);   % zero for angular velocity
R0_vgg = inv(K)*H_vgg*K;
angle0 = rotm2eul(R0_vgg./det(R0_vgg));
x0(1:3) = [angle0(3) angle0(2) angle0(1)];
x0(10) = f0;
%x0(1:3) = [0.1 -0.2 0.3];
% R0_bond = 3.14/2; w_bond = 0.2233/size(image1,1);
% lb = -1.*[R0_bond*ones(1,3), w_bond*ones(1,6).*v_scale];
% ub = 1.*[R0_bond*ones(1,3), w_bond*ones(1,6).*v_scale];
[x,resnorm,residual,exitflag,output,lambda,jacobian] = ...
lsqnonlin(@unk_f_rs_homo_cost_image,x0,lb,ub,options,A,b,K);

% H A1 A2 caculation 
disp('estimated H:')
H = R_caculate(x(1),x(2),x(3));
f0 = x(10);
K(1,1)=f0;
K(2,2)=f0;
H_image = K*H*inv(K);
scale = H_image(3,3);
H_image = H_image./scale;

w1_est = x(4:6);  A1_est = K*-H*X_(w1_est)*inv(K);  A1_est = A1_est./scale;
w2_est = x(7:9);  A2_est = K*X_(w2_est)*H*inv(K);   A2_est = A2_est./scale;

% RS image warping 
disp('RS warping begin...')
%generate warp map 
for row = 1 : size(image1, 1)
    for column = 1 : size(image1, 2)
        u1 = column;
        v1 = row;
        
        v2 = getv2(u1,v1,H_image,A1_est,A2_est);
        H_RS = H_image + A1_est*v1 + A2_est*v2;
        estimatedP = H_RS*[u1;v1;1];
        estimatedP = estimatedP(1:3)./estimatedP(3);
        uu(row,column) = estimatedP(1);
        vv(row,column) = estimatedP(2);
    end
end

% warp image range 
u_min = floor(min(min(uu))); 
u_max = ceil(max(max(uu))); 
v_min = floor(min(min(vv))); 
v_max = ceil(max(max(vv))); 

if u_min > 0
    u_min = 0;
end
if v_min > 0
    v_min = 0;
end

W = u_max - u_min;
H = v_max - v_min;
uu = uu - u_min;
vv = vv - v_min;
img1_warp = zeros(H,W,3);

% generate warp image 
image1 = im2double(image1);
for row = 1:size(image1, 1)
    for col = 1:size(image1, 2)
        u_warp = ceil(uu(row,col));
        v_warp = ceil(vv(row,col));
        u_warp_d = floor(uu(row,col));  
        v_warp_d = floor(vv(row,col));
        img1_warp(v_warp,u_warp ,:) = image1(row,col,:);
        if u_warp_d == 0
            u_warp_d = 1;
        end
        if v_warp_d == 0
            v_warp_d = 1;
        end
        img1_warp(v_warp_d,u_warp_d,:) = image1(row,col,:);
        img1_warp(v_warp_d,u_warp,:) = image1(row,col,:);
        img1_warp(v_warp,u_warp_d,:) = image1(row,col,:);
    end
end
% for i = 1:3
%     img1_warp(:,:,i) = imclose(img1_warp(:,:,i), strel('disk',1));
% end
figure, imshow(img1_warp)



% image merge 
% move B in respect to the top-left corner of A
mode = 'Average'
offsetW = round(-u_min);
offsetH = round(-v_min);
image_merge = blendMode(img1_warp, image2, mode, offsetW, offsetH);

figure(), imshow(image_merge); title('RS Stiching')
toc
% imgn = zeros(H,W,3);
% u_map = zeros(H,W);
% v_map = zeros(H,W);
% imgn(-v_min:(-v_min+size(image1,1)-1), -u_min:(-u_min+size(image1,2)-1),:) = double(image1);
% figure() ,imshow(imgn,[])
% u_map(-v_min:-v_min+size(image1,1)-1, -u_min:-u_min+size(image1,2)-1) = uu;
% v_map(-v_min:-v_min+size(image1,1)-1, -u_min:-u_min+size(image1,2)-1) = vv;
% 
% % generate RS warp image 
% [x,y] = ndgrid(0:(size(imgn,1)-1), 0:(size(imgn,2)-1)); % coordinate image
% % Interpolate updated image
% image1_warp = interpn(x,y,imgn(:,:,1),u_map,v_map,'linear',0); % moving image intensities at updated points
% figure() ,imshow(image1_warp',[])

%% image correction 
tic
disp('Stitching RS rectification begin...')
%generate warp map 
for row = 1 : size(image_merge, 1)
    for column = 1 : size(image_merge, 2)
        u = column ;
        v = row ;
        
        R0 = R_caculate(x(1),x(2),x(3));
        estimatedP = K*(eye(3) - X_(w1_est_inverse)*v)* inv(K) * [u;v;1];
        
        estimatedP = estimatedP./estimatedP(3);
        uu(row,column) = estimatedP(1);
        vv(row,column) = estimatedP(2);
    end
end

% warp image range 
u_min = floor(min(min(uu))); 
u_max = ceil(max(max(uu))); 
v_min = floor(min(min(vv))); 
v_max = ceil(max(max(vv))); 
W = u_max - u_min;
H = v_max - v_min;
uu = uu - u_min;
vv = vv - v_min;
img1_warp = zeros(H,W,3);

for row = 1:size(image_merge, 1)
    for col = 1:size(image_merge, 2)
        u_warp = ceil(uu(row,col));
        v_warp = ceil(vv(row,col));
        u_warp_d = floor(uu(row,col));  
        v_warp_d = floor(vv(row,col));
        img1_warp(v_warp,u_warp ,:) = image_merge(row,col,:);
        if u_warp_d == 0
            u_warp_d = 1;
        end
        if v_warp_d == 0
            v_warp_d = 1;
        end
        img1_warp(v_warp_d,u_warp_d,:) = image_merge(row,col,:);
        img1_warp(v_warp_d,u_warp,:) = image_merge(row,col,:);
        img1_warp(v_warp,u_warp_d,:) = image_merge(row,col,:);
    end
end
figure(), imshow(img1_warp); title('RS Stiching and Correction')
toc 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%       function we need                                                                                     
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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

v2 = real(v2_1);
if a == 0 
    v2 = 0;
end
end

function [u1 v1] = getu1v1(u2,v2,H,A1,A2)
end