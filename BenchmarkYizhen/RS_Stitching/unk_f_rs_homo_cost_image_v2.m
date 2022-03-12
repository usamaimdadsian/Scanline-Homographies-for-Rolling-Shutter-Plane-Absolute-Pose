function f = unk_f_rs_homo_cost_image_v2(x0,tracks,K,num,centerImageIdx, tforms)
f = [];

% prepare K matix
K(1,1) = x0(end);
K(2,2) = x0(end);

% prepare rotation and w
[R,w] = x0toParas(x0, num, centerImageIdx);

% counting errors
for i = 1:size(tracks,2)
    for j = 1:size(tracks(i).ViewIds,2)
        % q1 located in q1_viewId view and measured as q1
        q1_viewId = tracks(i).ViewIds(j);
        q1 = tracks(i).Points(j,:);
        for k = 1:size(tracks(i).ViewIds,2)
            if j ~= k
                % q2 located in q2_viewId view and measured as q2_measure
                q2_viewId = tracks(i).ViewIds(k);
                q2_measure = tracks(i).Points(k,:);
                % RS_H caculation
                R0 = R{q2_viewId}' * R{q1_viewId};
                H_GS = K * R0 * inv(K);
                A1 = K * -R0 * X_(w{q1_viewId}) * inv(K);
                A2 = K * X_(w{q2_viewId}) * R0 * inv(K);
                scale = H_GS(3,3);
                H_GS = H_GS./scale;
                A1 = A1./scale;
                A2 = A2./scale;
                
                % erros
                %q1_nor = inv(K) * [q1';1];
                q2_predic = get_q2(q1(1),q1(2),H_GS,A1,A2);
                %q2_predic = K * [q2_predic_nor;1];
                err1 = q2_predic(1) - q2_measure(1);
                err2 = q2_predic(2) - q2_measure(2);
                %             if abs(err1) > 5
                %                 err1 = 5;
                %             end
                %             if err2 > 5
                %                 err2 = 5;
                %             end
                f = [f; err1;err2];
            end
        end
    end
end
f = abs(double(f));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% RS homography translformation
function q2 = get_q2(u1,v1,H,A1,A2)
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

H_RS = H + A1*v1 + A2*v2;
estimatedP = H_RS*[u1;v1;1];
q2(1,1) = estimatedP(1)/estimatedP(3);
q2(2,1) = estimatedP(2)/estimatedP(3);
end

%% x0 to caemra parameters
function [R,w] = x0toParas(x0, num, centerImageIdx)
% rotation
R_vector = [x0(1:centerImageIdx*3-3) 0 0 0 x0(centerImageIdx*3-2:3*num)];
for i = 1:num
    R{i} = R_caculate(R_vector(i*3-2),R_vector(i*3-1),R_vector(i*3));
end
% angular velocity
w_vector = x0(3*num-2:end-1);
for i = 1:num
    w{i} = w_vector(i*3-2:i*3);
end
end