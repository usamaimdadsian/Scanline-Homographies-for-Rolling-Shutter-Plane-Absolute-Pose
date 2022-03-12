function [R] = R_caculate(x1,x2,x3)

R_1 = [1,0,0;
       0,cos(x1),-sin(x1);
       0,sin(x1),cos(x1)];
R_2 = [cos(x2),0,sin(x2);
       0,1,0;
       -sin(x2),0,cos(x2)];
R_3 = [cos(x3),-sin(x3),0;
       sin(x3),cos(x3),0;
       0,0,1];
R = R_3 * R_2 * R_1;  %rotation matrix 