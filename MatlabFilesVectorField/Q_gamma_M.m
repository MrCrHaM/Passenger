function Q = Q_gamma_M(Adj, u, gamma, r)
% Compute Q_gamma(M) = sum_{ij in E} M_ij L_ij - gamma sum_{i in V} M_ii L_ri
% where M \approx u*u' with n x 1 vector u
% Memory efficient, never keep non-sparse n x n matrices

n = size(Adj, 1);
degs = sum(Adj)';

%% method 1, requires computing n x n matrix u*u'
% tic; 
% X = u*u';
% Q = diag(sum(X.*Adj)) - X.*Adj;
% dX = degs.*diag(X);
% starX = diag(dX);
% starX(r,:) = -dX';
% starX(:,r) = -dX;
% starX(r,r) = sum(dX) - degs(r)*X(r,r);
% %     out = Q + ones(size(Q))*dX*ones(size(Q)) - gamma*diag(diag(X));
% Q1 = Q - gamma*starX;
% 
% toc
% 
% tic;

%rng(4);
%u = norminv(rand(n,1),0,1);

%% method 2, iterate over edges
[coords_i, coords_j] = ind2sub([n,n],find(triu(Adj)));
xvals = zeros(length(coords_i),1);
dvals = zeros(n,1);
for c = 1:length(coords_i);
    ii = coords_i(c);
    jj = coords_j(c);
    xvals(c) = -u(ii)*u(jj); %Y(coords_i(c),coords_i(c)) + Y(coords_j(c),coords_j(c)) - Y(coords_i(c),coords_j(c)) - Y(coords_j(c),coords_i(c));
    dvals(ii) = dvals(ii) + u(jj)/2;
    dvals(jj) = dvals(jj) + u(ii)/2;
end
dvals = dvals.*u;
Q = sparse([coords_i; (1:n)'], [coords_j; (1:n)'], [xvals; dvals], n, n);
Q = Q + Q';

dX = degs .* u.^2;
starX = diag(dX);
starX(r,:) = -dX';
starX(:,r) = -dX;
starX(r,r) = sum(dX) - degs(r)*u(r)^2;

Q = Q - gamma * starX;
% toc
% norm(Q1 - Q2, 'fro')