function out = paramQ_11(Adj, gamma, Y, r, X)
% implements connectivity constraint:
%   Q(M) = sum_{ij in E} M_ij L_ij + r 11^T - gamma diag(M)
% returns d/dX trace(Y*Q(X)) when called as paramQ_11(Adj, gamma, Y)
% returns Q(X)               when called as paramQ_11(Adj, gamma, [], r, X)

if nargin == 3
    % return d/dX trace(Y * Q(X))
    YY = zeros(size(Y));
    for i = 1:size(YY,1)
        for j = 1:i-1
            YY(i,j) = Adj(i,j) * (Y(i,i) + Y(j,j) - Y(i,j) - Y(j,i));
            YY(j,i) = YY(i,j);
        end
    end
    out = YY/2 - diag(diag(Y));
else
    % return Q(X)
    Q = zeros(size(X));
    for i = 1:size(Q,1)
        for j = i+1:size(Q,2)
            Lij = zeros(size(X));
            Lij(i,i) = 1; Lij(j,j) = 1; 
            Lij(i,j) = -1; Lij(j,i) = -1;
            Q = Q + Adj(i,j) * X(i,j) * Lij;
        end
    end
    out = Q + r*ones(size(Q)) - gamma*diag(diag(X));
%     out = Q + - gamma*diag(diag(X));
end