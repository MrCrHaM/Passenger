img = makeChessBoard(2,1,4);% makeChessBoard(20, 10, 4);
image(img);
colorbar;
A = makeGrid(img);
image(A)
sigma = 5;
e = 2;
graph = fully_connected_e_neighbour_graph(img, sigma, e);
[L, D] = laplacian(graph);
% inverse = inv(sparse(10 * L + eye(size(L, 1))));
proj = normrnd(0,1,[size(A, 1),1]);
c = (sparse(10 * L + eye(size(L, 1)))) \ proj;
gamma = 1;
s = zeros(1, size(L, 1))';
s(10) = 1;
s = 2;
k = 1;
max_iter = 500;
u_s = starOpt_fast(A, c, gamma, s, k, max_iter);