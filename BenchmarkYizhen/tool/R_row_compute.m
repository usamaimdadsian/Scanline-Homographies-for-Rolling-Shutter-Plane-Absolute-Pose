function R_row = R_row_compute(row, rowNum, R0_ato, a1, a2, b1, b2, c1, c2)

yf = (row-1)/rowNum;

rx = R0_ato(1) + a1*yf + a2*yf^2;
ry = R0_ato(2) + b1*yf + b2*yf^2;
rz = R0_ato(3) + c1*yf + c2*yf^2;

Z = 1 + rx^2 + ry^2 + rz^2;

R(1,1) = 1 + rx^2 - ry^2 - rz^2;
R(1,2) = 2*rx*ry - 2*rz;
R(1,3) = 2*ry + 2*rx*rz;
R(2,1) = 2*rz + 2*rx*ry;
R(2,2) = 1- rx^2 + ry^2 - rz^2;
R(2,3) = 2*ry*rz - 2*rx;
R(3,1) = 2*rx*rz - 2*ry;
R(3,2) = 2*rx + 2*ry*rz;
R(3,3) = 1 - rx^2 - ry^2 + rz^2;
R_row = (1/Z).*R;


