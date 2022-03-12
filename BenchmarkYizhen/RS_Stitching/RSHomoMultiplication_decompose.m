function H_RSn = RSHomoMultiplication_decompose(n, camParaPair, K)

R1 = camParaPair(n).R;
R0 = camParaPair(n-1).R;
R = R_caculate(R1(1), R1(2), R1(3)) * R_caculate(R0(1), R0(2), R0(3));
H_RSn.H_GS = K * R * inv(K);
H_RSn.A1 = K * (-R * X_(camParaPair(1).w)) * inv(K);
H_RSn.A2 = K *(X_(camParaPair(n).w) * R) * inv(K);

scale = H_RSn.H_GS(3,3);
H_RSn.H_GS  = H_RSn.H_GS ./ scale;
H_RSn.A1  = H_RSn.A1 ./ scale;
H_RSn.A2  = H_RSn.A2 ./ scale;
end