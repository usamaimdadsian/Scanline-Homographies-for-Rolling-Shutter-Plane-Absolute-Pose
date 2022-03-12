function f = unk_f_rs_homo_cost_image(x0,A,b,K)

f0 = x0(10);
K(1,1)=f0;
K(2,2)=f0;

R0 = R_caculate(x0(1),x0(2),x0(3));
w1 = x0(4:6);  A1 = K*-R0*X_(w1)*inv(K);
w2 = x0(7:9);  A2 = K*X_(w2)*R0*inv(K);
R0 = K * R0 * inv(K);
scale = R0(3,3);
R0 = R0./scale;
A1 = A1./scale;
A2 = A2./scale;

x1 = R0'; x2 = A1'; x3 = A2';
x = [x1(:);x2(:);x3(:)]; 
x(9) = [];

f = A*x - b;
%f(end+1) = R0(3,3) -1;
f = double(f);
