clear all;
close all;
cvx_solver sedumi

%% data generation, generate n1xn2 lattice, 4-connected
% n1 = 5; n2 = 5;  % smaller example
n1 = 8; n2 = 8;  % larger example
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
% xx_line = zeros(n,1);
% xx_line([7, 9, 10, 15]) = 1;  % 4 nodes with value 1, 1 node disconnected
% xx = reshape(xx_line, n1, n2);

% larger example
xx(3:4, 2:4) = 1;  % larger connected component
xx(5, 6:8) = 1;  % smaller connected component
K = 9;

% observations with some noise
yy = xx + 0.01*randn(size(xx));

% plot using imagesc to display as matrix
figure, imagesc(xx), title('Original');

%% cvx primal opt
M = subgraphOpt(Adj, yy(:), K, 0.1);

% show result using imagesc, no thresholding
figure, imagesc(reshape(diag(M), n1, n2)), colorbar, title('Estimated');

%% project and display
S = diag(M) > 1e-2;
Msub = M(find(S), find(S));  % non-zero submatrix of M

% lattice coordinates of indices in M
[ind1, ind2] = ind2sub([n1 n2], (1:n)');
ind_text = strcat(num2str(ind1), ',', num2str(ind2));
ind_sub = ind_text(S, :);

% random projection to 3-dim
% R = (rand(size(Msub,1),3) > 0.5)*2 - 1;  % random binary vectors
R = randn(size(Msub,1),3);  % random Gaussian vectors
R = R./repmat(sqrt(sum(R.^2)), size(Msub,1), 1);  % normalize to unit length
V = Msub*R;
figure, scatter3(V(:,1), V(:,2), V(:,3)), hold on
text(V(:,1), V(:,2), V(:,3), ind_sub), title('Random projection')
