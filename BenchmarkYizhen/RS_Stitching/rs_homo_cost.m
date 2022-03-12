function f = rs_homo_cost(x0,A,b)

R0 = R_caculate(x0(1),x0(2),x0(3));
R0 = R0 ./ R0(3,3);
w1 = x0(4:6);  A1 = -R0*X_(w1);
w2 = x0(7:9);  A2 = X_(w2)*R0;

x1 = R0'; x2 = A1'; x3 = A2';
x = [x1(:);x2(:);x3(:)]; 
x(9) = [];

f = A*x - b;
%f(end+1) = R0(3,3) -1;
f = double(f);
