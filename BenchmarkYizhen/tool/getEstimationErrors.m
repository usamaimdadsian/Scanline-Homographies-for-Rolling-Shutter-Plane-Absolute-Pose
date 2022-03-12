function [R_error,t_error,w_error,d_error] = getEstimationErrors(R_est,t_est,w_est,d_est,R_GT,t_GT,w_GT,d_GT)

R_error = angle2vecotrs(R_est*[0 0 1]',R_GT*[0 0 1]');
t_error = norm(t_est - t_GT);
w_error = norm(w_est - w_GT);
d_error = norm(d_est - d_GT);

function angle = angle2vecotrs(u,v)
angle = atan2d(norm(cross(u',v')),dot(u',v'));
