function [img1_warp, u_min, v_min] = H_RS_ImageWarping(image1, H_image, A1_est, A2_est, K, H_RS_inverse)

%generate warp map 
for row = 1 : size(image1, 1)
    for column = 1 : size(image1, 2)        
        H_RS = H_image + A1_est*row + A2_est*getv2(column,row,H_image,A1_est,A2_est);
        estimatedP = H_RS*[column;row;1];
        uu(row,column) = estimatedP(1)/estimatedP(3);
        vv(row,column) = estimatedP(2)/estimatedP(3);
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
img1_warp = zeros(H,W,3);
% img1_warp = ones(H,W,3);

uu = uu - u_min;
vv = vv - v_min;
uu = ceil(uu);
vv = ceil(vv);

% generate warp image 
% image1 = im2double(image1);
% for row = 1:size(image1, 1)
%     for col = 1:size(image1, 2)
%         u_warp = uu(row,col);
%         v_warp = vv(row,col);
%         img1_warp(vv(row,col),uu(row,col) ,:) = image1(row,col,:);
%         %img1_warp(vv(row,col)+1,uu(row,col) + 1 ,:) = image1(row,col,:);
%     end
% end

% dense mapping 
H_GS = H_RS_inverse.H_GS;
A1 = H_RS_inverse.A1;
A2 = H_RS_inverse.A2;

image1 = im2double(image1);
[row_max, col_max] = size(image1(:,:,1));
for row = v_min+1:size(img1_warp, 1)
    for col = u_min+1:size(img1_warp, 2)
        % use [u2 v2] to get [u1 v1] 
        % it is a inverse mapping 
        v1 = getv2(col,row,H_GS,A1,A2);
        p = (H_GS + A1*row + A2*v1)*[col;row;1];
        u1 = p(1)/p(3);
        
        if (v1>0 && v1<row_max && u1>0 && u1<col_max)
            % img1_warp(row - v_min,col - u_min,:) = image1(ceil(v1),ceil(u1),:); 
            img1_warp(row - v_min,col - u_min,:) = reSample(u1, v1, image1); 
        end
    end
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function v2 = getv2(u1,v1,H,A1,A2)
p = [u1;v1;1];
a = A2(3,:)*p;
b = H(3,:)*p + A1(3,:)*p*v1 - A2(2,:)*p;
c = -(H(2,:)+A1(2,:)*v1)*p;

v2_1 = (-b + sqrt(b^2 - 4*a*c))/(2*a);
v2_2 = (-b - sqrt(b^2 - 4*a*c))/(2*a);

med = 500;
dis1 = abs(v2_1 - med);
dis2 = abs(v2_2 - med);
if dis1 < dis2 
    v2 = v2_1;
else
    v2 = v2_2;
end
v2 = real(v2);
if a == 0 
    v2 = 0;
end
end

%% re-sampling of image pixel 
function pixelValue = reSample(u, v, image)
u_c = round(u);
v_c = round(v);

if (u_c > 0 && v_c>0)
    pixelValue = image(v_c,u_c,:);
else
    pixelValue = image(ceil(v),ceil(u),:);
end
end