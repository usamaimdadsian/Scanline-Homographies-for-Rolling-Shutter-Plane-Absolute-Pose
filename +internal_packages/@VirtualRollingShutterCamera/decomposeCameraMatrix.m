% [K, R, C] = decomposeCameraMatrix (P)
% A function to decomposition camera matrix P to K * [R, -R*C]
%    - K   Camera Instrinsics 
%    - R   Camera Rotation
%    - C   Camera Center Position

function [K, R, C] = decomposeCameraMatrix (P)

[~, ~, V] = svd(P);
hgc = V(:, end);

C = hgc([1,2,3]) ./ hgc(4);

M = P(:, [1,2,3]);

[K, R] = RQ_Decomposition (M);

% the diagonal elements of K are positive
for ii = [1, 2, 3]
    if (K(ii, ii) < 0)
        K(:, ii) = - K(:, ii);
        R(ii, :) = - R(ii, :);
    end
end

if ( norm(K*R * [eye(3), -C] - P, 'fro')  > 1e-8)
    fprintf(2, '\nLarge error in camera matrix decomposition: %f\n', norm(K*R * [eye(3), -C] - P, 'fro'));
    fprintf(2, 'det(P33) = %.18f \n', det(M));
%     fprintf(2, 'P(3,:) = [%f   %f   %f   %f] \n', P(3,:));
%     fprintf(2, 'Error in RQ.  norm(K*R - P33) = %.18f \n', norm(K*R-M));
%     fprintf(2, 'Error in P * hgc.  norm(P*hgc) = %1.18f \n', norm(P*hgc));
%     fprintf(2, 'hgc = [%f   %f   %f   %.18f] \n', hgc(1), hgc(2), hgc(3), hgc(4));
%     fprintf(2, 'Error in P [C; 1].  norm(P*[C;1]) = %.18f \n', norm(P*[C;1]));
%     fprintf(2, '\n');
end

if(det(R) < 0)
    R = -R;
end

%K = K/K(3,3);

end



% RQ decomposition from classical QR decomposition. A must be a square matrix
% A = R * Q.   R - upper traingular matrix,  Q - orthogonal matrix
function [R, Q] = RQ_Decomposition (A)

if (size(A,1) ~= size(A,2))
    fprintf(2, 'Error: @RQ_Decomposition, the input matrix must be square\n');
end

nsize = length(A);

P = eye(nsize);

P = P(:,  end:-1:1);

PA = P * A;

[Q, R] = qr(PA');

Q = P * Q';

R = P * R' * P;

end
