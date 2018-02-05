function fn = jing_SDP(Adj, params, c, K, gamma, p)

% Performs maximization over subgraphs of graph given by its adjacency matrix Adj with node costs given by c.
%   Adj: adjacency matrix of input graph with n nodes, n x n
%   c: node costs to be maximized, n x 1
%   K: sparsity constraint on chosen nodes, none if K = 0
%   gamma: internal conductance of subgraph to be found
%   M: n x n matrix solution

% C = c*c';
n = length(c);

% cvx_begin
%     variable M(n,n) symmetric  % main variable
%     expression L_M(n,n)  % laplacian of M
%     expression A_M(n,n)  % A.M
%     expression L_A_M(n,n)  % laplacian of A.M
% 
%     L_M = diag(sum(M,2)) - M;
%     A_M = Adj.*M;
%     L_A_M = diag(sum(A_M,2)) - A_M;
% 
%     maximize( c'*diag(M) )
%     subject to
%         M == semidefinite(n)
%         M >= 0
%         trace( M ) <= K
%         L_A_M - gamma*L_M == semidefinite(n)
% 
%         % constaints from NIPS paper
%         M <= 1
%         M(p,p) == 1
%         diag(M) <= M(:,p)
%         for i = 1:n
%             M(i,:) <= M(i,i)
%         end
% cvx_end

cvx_begin 
    variables fn(n) fe(params.m);    % node / edge indicator variables

    minimize( -c' * fn ); % + 0.3*( sum(fe) );     % regularization on edge sum
    subject to
       params.DI_AoM' * diag([fn; fe]) * params.DI_AoM - gamma * params.DI_M' * diag([fn; fe]) * params.DI_M == semidefinite(n);       
       fn >= 0;
       fn <= 1;
       fe >= 0;
       fe <= params.DI_p1*fn;
       fe <= params.DI_p2*fn;
       fn(p) == 1;       
       sum(fn) <= K;        %  size constraint
cvx_end

