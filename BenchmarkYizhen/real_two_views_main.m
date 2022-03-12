% Real RS images stitching  
% Author: Yizhen LAO
% Two views example 
addpath('./data');
addpath('./RS_Stitching');
addpath(genpath('./tool'));

clc 
close all
clear all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% camera calibration info
% rows=480; 
% cols=640;
% f = 35 * cols / (25.4 * 0.98);
% K = [f 0 cols/2;
%     0 f rows/2;
%     0 0 1];

% Univ. Adelaide Camera Calib.
% fx = 481.203;
% fy = 480.0;
% f = (fx+fy)/2;
% cx = 320.5;
% cy = 240.5;
% K = [fx 0 cx;
%     0 fy cy;
%     0 0 1];

% LiU real Camera Calib.
% fx = 1505.1283359786307;
% fy = 1513.7789208311444;
% f = (fx+fy)/2;
% cx = 657.81734686405991;
% cy = 349.91807538147589;

% Blender synthetic camera
rows=1080; 
cols=1920;
fx = 1.9177e+03;
fy = 1.9034e+03;
f = (fx+fy)/2;
cx =960.0;
cy =540.0;



% K = [fx 0 cx;
%     0 fy cy;
%     0 0 1];


% % % % data reading 
% % % path = './data/159/';
% image1 = imread('E:\dataset\RollingShutter\UnivAdelaide_dataset\Sequence77\Sequence77\RS_image\RS_0001.png');
% image2 = imread('E:\dataset\RollingShutter\UnivAdelaide_dataset\Sequence77\Sequence77\RS_image\RS_0036.png');
% figure(), imshow(image1); title('Original input 1')
% figure(), imshow(image2); title('Original input 2')
% op_img_path = strcat('E:\RESULTS\RollingShutter_CVPR2021\UnivAdelaide\Yizhen2\');

% data reading 
path = './data/159/';
image1 = imread('E:\dataset\RollingShutter\UCA\Real\SEQ1\frame132_m.png');
image2 = imread('E:\dataset\RollingShutter\UCA\Real\SEQ1\frame1037_m.png');
figure(), imshow(image1); title('Original input 1')
figure(), imshow(image2); title('Original input 2')
op_img_path = strcat('E:\RESULTS\RollingShutter_CVPR2021\UCA_Book\Frame1037\Yizhen2View\');

gs_warped_on_rs_path = strcat(op_img_path,'gs2rs.png');
rs_warped_on_gs_path = strcat(op_img_path,'rs2gs.png');
rs_rectified_path = strcat(op_img_path,'rsRectified.png');


%%% This is for iPad Air 2
focalLength    = [1.9177e+03 1.9034e+03]; % focal length (fx,fy) of rolling-shutter camera
principalPoint =[956.2334 537.0888]; % principal point of rolling-shutter camera
radialDistortion = [-0.0154 0.9590]; % radial distortion of rolling-shutter camera
imageSize  = [1080 1920]; % image height, width of rolling-shutter image
intrinsics  = [focalLength(1) 0 0; 0 focalLength(2) 0; principalPoint(1), principalPoint(2) 1];
intrinsics = cameraParameters('IntrinsicMatrix',intrinsics,'RadialDistortion',radialDistortion);
image1 = undistortImage(image1,intrinsics);
image2 = undistortImage(image2,intrinsics);

K = [focalLength(1) 0 principalPoint(1);
    0 focalLength(2) principalPoint(2);
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
figure; ax = axes;
showMatchedFeatures(I1,I2,matchedPoints1,matchedPoints2,'montage','Parent',ax);

save('DataReal1_iPad_book.mat','image1','image2','IP1','IP2','K');

%% GS image transformation 
disp('[GS-based Homography...............................................................')
fprintf('Image1 size (%d, %d)\n',size(image1,1),size(image1,2));
[H_vgg, inliers] = ransacfithomography_vgg([IP1;ones(1,matchedPoints1.Count)],[IP2;ones(1,matchedPoints2.Count)], 0.01);
size(image1)
size(H_vgg)
H_vgg
inliers
[image1_warp, R1] = imwarp(image1,projective2d(H_vgg'));
size(image1_warp)
R2 = imref2d(size(image2));
%image_merge = imfuse(image1_warp, R1, image2, R2, 'method', 'blend');
mode = 'Average';
figure(), imshow(image1_warp); title('Pre GS Stitching');
pause(0.1);
% move B in respect to the top-left corner of A
offsetW = round(-R1.XWorldLimits(1));
offsetH = round(-R1.YWorldLimits(1));
image_merge = blendMode(image1_warp, image2, mode, offsetW, offsetH);

figure(), imshow(image_merge); title('GS Stitching')
disp('********************************************************************************************************************')
imwrite(image_merge,gs_warped_on_rs_path);
pause(0.01);
close all;
%% RS image transformation 

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
x0 = zeros(1,9);   % zero for angular velocity
R0_vgg = inv(K)*H_vgg*K;
angle0 = rotm2eul(R0_vgg./det(R0_vgg));
x0(1:3) = [angle0(3) angle0(2) angle0(1)];
%x0(1:3) = [0.1 -0.2 0.3];
R0_bond = 3.14/2; w_bond = 0.2233/size(image1,1); v_scale = [0.2 0.2 1 0.2 0.2 1];
lb = -1.*[R0_bond*ones(1,3), w_bond*ones(1,6).*v_scale];
ub = 1.*[R0_bond*ones(1,3), w_bond*ones(1,6).*v_scale];
[x,resnorm,residual,exitflag,output,lambda,jacobian] = ...
lsqnonlin(@rs_homo_cost_image,x0,lb,ub,options,A,b,K);

% H A1 A2 caculation 
disp('estimated H:')
H = R_caculate(x(1),x(2),x(3));
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
x0 = zeros(1,9);   % zero for angular velocity
R0_vgg = inv(K)*H_vgg*K;
angle0 = rotm2eul(R0_vgg./det(R0_vgg));
x0(1:3) = [angle0(3) angle0(2) angle0(1)];
%x0(1:3) = [0.1 -0.2 0.3];
R0_bond = 3.14/2; w_bond = 0.2233/size(image1,1);
lb = -1.*[R0_bond*ones(1,3), w_bond*ones(1,6).*v_scale];
ub = 1.*[R0_bond*ones(1,3), w_bond*ones(1,6).*v_scale];
[x,resnorm,residual,exitflag,output,lambda,jacobian] = ...
lsqnonlin(@rs_homo_cost_image,x0,lb,ub,options,A,b,K);

% H A1 A2 caculation 
disp('estimated H:')
H = R_caculate(x(1),x(2),x(3));
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
figure, imshow(img1_warp); title('RS after warping')
imwrite(img1_warp,rs_rectified_path);
% ^^^ This is the rectified image (-Agniva)

% image merge 
% move B in respect to the top-left corner of A
mode = 'Average'
offsetW = round(-u_min);
offsetH = round(-v_min);
image_merge = blendMode(img1_warp, image2, mode, offsetW, offsetH);

figure(), imshow(image_merge); title('RS Stiching')


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


imwrite(img1_warp,rs_warped_on_gs_path);

close all;

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