function [A,b] = GenerateCroutMatrix(p)% p is the size of the inner Aij matrices, the matrix A will be of size 3*p by 3*p
%Generates a random Crout matrix that has the proper form according to
%smoks_help.jpg and the_question.jpg 
a12 = rand(p);
a23 = rand(p);
A = [eye(p)*3*p,a12,zeros(p); transpose(a12),eye(p)*3*p, a23;zeros(p),transpose(a23),eye(p)*3*p];
b = transpose(1:3*p);
end

