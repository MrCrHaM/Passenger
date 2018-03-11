img = makeChessBoard(20, 10, 3);% makeChessBoard(20, 10, 4);
imagesc(img);
colorbar;
A = makeGrid(img);

sigma = 5;
e = 2;
graph = fully_connected_e_neighbour_graph(img, sigma, e);
[L, D] = laplacian(graph);
max_iter = 30;
% inverse = inv(sparse(10 * L + eye(size(L, 1))));
rng(1);
% proj1 = randn([size(A,1),1]);
d = 15;
proj1 = norminv(rand(size(A,1),d),0,1);
c = (sparse(10 * L + eye(size(L, 1)))) \ proj1;
gamma = .01;
rowN = round(size(img,1) / 2);
colN = round(size(img,2) / 2);
s = rowN * (size(img,2)) + colN;
k = 20;

u_s = starOpt_fast(A, inv(D) * c, gamma, s, k, max_iter);
% proj2 = randn([max_iter,1]);
rng(3 + max_iter);
proj2 = norminv(rand(max_iter,1),0,100);

v=diag(u_s*u_s');
%v=diag(u_s*proj2*proj2'*u_s');
r = reshape(v, sqrt(size(v,1)),sqrt(size(v,1)));
r(rowN, colN) = r(rowN, colN)/2;
disp(max(max(r)))
if(rowN ~= 1) 
    r(rowN - 1, colN) = r(rowN - 1, colN)/2; 
end
disp(max(max(r)))
if(rowN ~= size(img,1)) 
    r(rowN + 1, colN) = r(rowN + 1, colN)/2; 
end
disp(max(max(r)))
if(colN ~= 1) 
    r(rowN, colN - 1) = r(rowN, colN - 1)/2; 
end
disp(max(max(r)))
if(colN ~= size(img,2)) 
    r(rowN, colN + 1) = r(rowN, colN + 1)/2; 
end
disp(max(max(r)))
imagesc(r); colormap('parula'); colorbar;
% imagesc(r); colormap('hot'); colorbar;

