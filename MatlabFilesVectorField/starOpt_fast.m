function u_s = starOpt_fast(A, c, gamma, s, k, max_iter)

n = size(A, 1);
A = sparse(A);

[coords_i, coords_j] = ind2sub([n, n], find(A));

[tri_i, tri_j] = ind2sub([n, n], find(triu(A)));
map = triuToFullIdx(coords_i, coords_j, tri_i, tri_j, A);

degs = sum(A)';
% degs = ones(n, 1);

eta = 5;
p = 2*(c'*c)/gamma;

tic
% randomly init Y ~~ y y' and normalize weighted trace
y = randn(n, k);
y = y * sqrt( p / ( sum(sum( degs(:, ones(k,1)) .* (y.^2) )) ) );

u_s = zeros(n, max_iter);
grad = sparse(n, n);

opts.issym = 1;
opts.tol = 1/n;
opts.p = 5;
opts.maxit = 100;
for i = 1:max_iter
    % define operator mtxop(x) = A*x for eigs(A)
    % where A = c*c' + P_gamma(Y)
    ydist = precomputeDists(tri_i, tri_j, map, y);
    mtxop = @(x) cc_P_gamma_Y_mult(coords_i, coords_j, degs, ...
                                   c, y, gamma, s, x, ydist);

    [u, ~] = eigs(mtxop, n, 1, 'la', opts);
    opts.v0 = u;
    u_s(:,i) = u;

    QX = Q_gamma_M(A, u, gamma, s);

    % MWU for dual variable Y
    grad = grad - QX;
    y = expleja(eta/2, grad, randn(n, k)/sqrt(k), [0, 1/n]);

    % normalize such that y*y' has degree weighted trace p
    y = y * sqrt( p / ( sum(sum( degs(:, ones(k, 1)) .* (y.^2) )) ) );
end
toc
