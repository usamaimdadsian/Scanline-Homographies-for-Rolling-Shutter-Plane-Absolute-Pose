
function [x] = SolveByCrout(A,b)
[m,n] = size(A);
if m ~= n
    disp('A is not an n x n matrix! (Is not square)')
    return
end
if rem(n,3) ~= 0
     disp('The property n = 3*p where n and p are natural numbers does not hold therefore A is an invalid input matrix')
    return
end
p = n / 3;
[L,U] = CroutPartition(A,p);
%LUx = b...Solving Ax = b in matlab is done by x = A\b; NOT by x = b/A
y = L\b;%first we solve Ly = b
%then we can substitute y into Ux = y
x = U\y;


Xex = A \ b % The exact solution?
disp('The relative error of the computed solution in percent');
str = norm(x - Xex)/norm(Xex)*100;% by default the norm function gives the 2-norm or eucledian norm Source: http://www.mathworks.com/help/matlab/ref/norm.html
disp(str);

end


