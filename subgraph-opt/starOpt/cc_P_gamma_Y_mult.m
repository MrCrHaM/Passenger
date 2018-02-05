function v = cc_P_gamma_Y_mult(coords_i, coords_j, degs, c, y, gamma, r, x, ydist)
% Compute (c*c' + P_gamma(Y)) * x, where 
% P_gamma(Y) = sum_{i,j in E} (L_ij * Y) e_ij - gamma sum_i d_i (L_ri * Y) e_ii
% and Y \approx y*y' with n x k matrix y
% ydist(e) = ||y(i,:) - y(j,:)||^2 for e = (i,j) in 1,...,m

n = length(degs);
assert(length(x) == n && size(y,1) == n);

% m = length(coords_i);
% 
% v = (c'*x) * c;
% for ind = 1:m
%     ii = coords_i(ind);
%     jj = coords_j(ind);
% 
%     v(ii) = v(ii) + x(jj) * ydist(ind);
% %     v(jj) = v(jj) + x(ii) * ydist(ind);
% %     v([ii jj]) = v([ii jj]) + [x(jj); x(ii)] * sum((y(ii,:) - y(jj,:)).^2);
% %     v([ii jj]) = v([ii jj]) + [x(jj); x(ii)] * ydist(ind);
% end
% v = v - gamma * (degs .* x .* sum((repmat(y(r,:), n, 1) - y).^2, 2));

zz = ydist .* x(coords_j);
vv = accumarray(coords_i, zz, [n 1]);
v = (c'*x) * c + vv - gamma * (degs .* x .* sum((repmat(y(r,:), n, 1) - y).^2, 2));
