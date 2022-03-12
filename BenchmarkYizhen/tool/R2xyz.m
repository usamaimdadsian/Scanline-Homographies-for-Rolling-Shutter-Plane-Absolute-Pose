function r = R2xyz(R)

z = (R(2,1) - R(1,2))/2;
y = (R(1,3) - R(3,1))/2;
x = (R(3,2) - R(2,3))/2;

r = [x y z];