clear all;
close all;
cvx_solver sedumi

%% data generation, generate n1xn2 lattice, 4-connected
n1 = 5; n2 = 5;  % smaller example
% n1 = 6; n2 = 8;  % larger example
n = n1*n2;

% compute unweighted adjacency matrix, Laplacian matrix
Adj = zeros(n, n);
for i=1:n1
    for j=1:n2
        if j<n2
            Adj((i-1)*n2+j, (i-1)*n2+j+1) = 1;
            Adj((i-1)*n2+j+1, (i-1)*n2+j) = 1;
        end
        if i<n1
            Adj((i-1)*n2+j, i*n2+j) = 1;
            Adj(i*n2+j, (i-1)*n2+j) = 1;
        end
    end
end

xx = zeros(n1, n2);

% smaller example
xx_line = zeros(n,1);
xx_line([7, 9, 10, 15]) = 1;  % 4 nodes with value 1, 1 node disconnected
xx = reshape(xx_line, n1, n2);

% larger example
% xx(3:4, 2:4) = 1;  % larger connected component
% xx(5, 6:8) = 1;  % smaller connected component

% observations with some noise
yy = xx + 0.01*randn(size(xx));
yy(yy < 0) = 0;

% plot using imagesc to display as matrix
% figure, imagesc(xx), title('Original');

% A_full = ones(n)-eye(n);
% L_KV = diag(sum(A_full)) - A_full;

%% parameters
gamma = 0.5;

c = yy(:);
C = c*c';
r = 1e-2;
p = 2*norm(C,'fro')/gamma;

%% cvx solution
tic
A = Adj;

cvx_begin
    variable M(n,n) symmetric  % main variable
    variable s
    expression L_M(n,n)  % laplacian of M
    expression A_M(n,n)  % A.M
    expression L_A_M(n,n)  % laplacian of A.M

    A_M = Adj.*M;
    L_A_M = diag(sum(A_M,2)) - A_M;

    maximize( trace(C*M) - p*s )
    subject to
        M == semidefinite(n)
        M >= 0
        trace( M ) == 1
        L_A_M + r*ones(n) - gamma*diag(diag(M)) + s*eye(n) == semidefinite(n)
cvx_end
toc

% show result using imagesc, no thresholding
figure, imagesc(reshape(diag(M), n1, n2)), colorbar, title('Estimated');

%% MWU solution
max_iter = 1000;
eta = 0.05;

tic
% randomly init Y
Y = (eye(n)+0.1*rand(n)); Y = Y'*Y; Y = Y/trace(Y)*p;
X = zeros(n,n,max_iter);

for i=1:max_iter
    % first eigenvector matrix for primal opt
    [u,v] = eigs(real(paramQ_11(Adj,gamma,Y) + C), 1, 'la');
    X(:,:,i) = u*u';
    QX = paramQ_11(Adj,gamma,[],r,X(:,:,i));
    
    % MWU for dual variable Y
    Y = expm(logm(Y) - eta*QX);
    
    Y = Y/trace(Y)*p;  % normalize to trace p
end
final = mean(X,3);
toc

% Qf = paramQ_11(Adj,gamma,[],r,final);
% ee = min(eig(Qf));

figure, imagesc(reshape(diag(final), n1, n2)), colorbar, title('Estimated');
% figure, imagesc(final), colorbar
