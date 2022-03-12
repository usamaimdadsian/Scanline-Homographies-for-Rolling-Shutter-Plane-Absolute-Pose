function H_RS_13 = RSHomoMultiplication(H_RS_12, H_RS_23)

H_RS_13.H_GS = H_RS_23.H_GS * H_RS_12.H_GS;
H_RS_13.A1 = H_RS_23.H_GS * H_RS_12.A1;
H_RS_13.A2 = H_RS_23.A2 * H_RS_12.H_GS;

scale = H_RS_13.H_GS(3,3);
H_RS_13.H_GS  = H_RS_13.H_GS ./ scale;
H_RS_13.A1  = H_RS_13.A1 ./ scale;
H_RS_13.A2  = H_RS_13.A2 ./ scale;
end