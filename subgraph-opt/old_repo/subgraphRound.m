function subgraphRound(M, n1, n2)

n = n1*n2;
% show result using imagesc, no thresholding
figure, imagesc(reshape(diag(M), n1, n2)), colorbar, title('Estimated');

% project and display
S = diag(M) > 1e-2;
Msub = M(find(S), find(S));  % non-zero submatrix of M

% lattice coordinates of indices in M
[ind1, ind2] = ind2sub([n1 n2], (1:n)');
ind_text = strcat(num2str(ind1), ',', num2str(ind2));
ind_sub = ind_text(S, :);

% random projection to 2-dim
% R = (rand(size(Msub,1),2) > 0.5)*2 - 1;  % random binary vectors
R = randn(size(Msub,1),2);  % random Gaussian vectors
R = R./repmat(sqrt(sum(R.^2)), size(Msub,1), 1);  % normalize to unit length
V = Msub*R;
figure, scatter(V(:,1), V(:,2)), hold on
text(V(:,1), V(:,2), ind_sub), title('Random projection')
