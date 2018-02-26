function out = paramQ_star(Adj, degs, gamma, Y, r, X)
% implements connectivity constraint:
%   Q(M) = sum_{ij in E} M_ij L_ij - gamma sum_{i in V} M_ii L_ri
% returns d/dX trace(Y*Q(X)) when called as paramQ_star(Adj, degs, gamma, Y, r)
% returns Q(X)               when called as paramQ_star(Adj, degs, gamma, [], r, X)

if nargin == 5
    % return d/dX trace(Y * Q(X))
    n = size(Y,1);
    [coords_i, coords_j] = ind2sub([n,n],find(triu(Adj)));
    yvals = zeros(length(coords_i),1);
    for c = 1:length(coords_i);
        yvals(c) = Y(coords_i(c),coords_i(c)) + Y(coords_j(c),coords_j(c)) - Y(coords_i(c),coords_j(c)) - Y(coords_j(c),coords_i(c));
    end
    YY = sparse(coords_i,coords_j,yvals,n,n);
    YY = YY + YY';
%     out = YY - gamma*spdiags(Y(r,r) + diag(Y) - Y(r,:)' - Y(:,r), 0, n, n);
    out = YY - gamma*spdiags(degs.*(Y(r,r) + diag(Y) - Y(r,:)' - Y(:,r)), 0, n, n);
else
    % return Q(X)
    Q = diag(sum(X.*Adj)) - X.*Adj;
    dX = degs.*diag(X);
    starX = diag(dX);
    starX(r,:) = -dX';
    starX(:,r) = -dX;
    starX(r,r) = sum(dX) - degs(r)*X(r,r);
%     out = Q + ones(size(Q))*dX*ones(size(Q)) - gamma*diag(diag(X));
    out = Q - gamma*starX;
end