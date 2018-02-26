img = makeChessBoard(20, 10, 4);
image(img);
colorbar;
A = makeGrid(img);
sigma = 5;
e = 2;
graph = fully_connected_e_neighbour_graph(img, sigma, e);
[L, D] = laplacian(graph);
max_iter = 500;
% inverse = inv(sparse(10 * L + eye(size(L, 1))));
proj = normrnd(0,1,[size(A, 1),1]);
c = (sparse(10 * L + eye(size(L, 1)))) \ proj;
gamma = 1;
s = zeros(1, size(L, 1))';
s(10) = 1;
s = 2000;
k = 1;

u_s = starOpt_fast(A, inv(D)* c, gamma, s, k, max_iter);

v = u_s*proj;
r = reshape(v, size(v, 1) ^ .5, size(v, 1) ^ .5);
image(img2)
hold on
image(img)
image(r)

proj = normrnd(0,1,[max_iter,1]);

v=u_s*proj;
r=reshape(v, sqrt(size(v,1)),sqrt(size(v,1)));
image(1000 * r)
