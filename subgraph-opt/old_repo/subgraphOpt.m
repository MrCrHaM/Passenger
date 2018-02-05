function M = subgraphOpt(Adj, c, K, gamma)

% Performs maximization over subgraphs of graph given by its adjacency matrix Adj with node costs given by c.
%   Adj: adjacency matrix of input graph with n nodes, n x n
%   c: node costs to be maximized, n x 1
%   K: sparsity constraint on chosen nodes, none if K = 0
%   gamma: internal conductance of subgraph to be found
%   M: n x n matrix solution

C = c*c';
n = length(c);

cvx_begin
    variable M(n,n) symmetric  % main variable
    expression L_M(n,n)  % laplacian of M
    expression A_M(n,n)  % A.M
    expression L_A_M(n,n)  % laplacian of A.M

    L_M = diag(sum(M,2)) - M;
    A_M = Adj.*M;
    L_A_M = diag(sum(A_M,2)) - A_M;

    maximize( trace( C * M ) )
    subject to
        M == semidefinite(n)
        M >= 0
        trace( M ) == 1
        (K > 0)*sum(sum( M )) <= K
        L_A_M - gamma*L_M == semidefinite(n)

        % constaints from NIPS paper
%         M <= 1
%         M(p,p) == 1
%         for i = 1:n
%             M(i,i) <= M(p,i)
%             M(i,:) <= M(i,i)
%         end
cvx_end

