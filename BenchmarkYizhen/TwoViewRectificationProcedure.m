function [rectification] = TwoViewRectificationProcedure (IP1, IP2, image1, image2, K)

%% GS image transformation 
disp('[GS-based Homography...............................................................')
[H_vgg, inliers] = ransacfithomography_vgg([IP1;ones(1,size(IP1,2))],[IP2;ones(1,size(IP2,2))], 0.01);

  %% RS image transformation 

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
            H_landmrk{row,column} = H_RS;
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
   % figure, imshow(img1_warp); title('RS after warping');
    
    rectification = img1_warp;
    % ^^^ This is the rectified image (-Agniva)


end

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
        
        if v2_1 < 0 || v2_1 >2000
            v2 = v2_2;
        else
            v2 = v2_1;
        end
        
        v2 = real(v2_1);
        if a == 0
            v2 = 0;
        end
    end