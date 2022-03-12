function H_RS_21 = RSHomoInverse(H_RS_12, K)

H_RS_21.H_GS = K * (inv(K) * H_RS_12.H_GS * K)' * inv(K);
H_RS_21.A1 = K * (inv(K) * H_RS_12.A2 * K)' * inv(K);
H_RS_21.A2 = K * (inv(K) * H_RS_12.A1 * K)' * inv(K);
