function map = triuToFullIdx(coords_i, coords_j, tri_i, tri_j, A)
% Finds a m x 1 mapping such that for any edge e = (i,j) in the full Adj
% matrix, map(e) returns the corresponding edge coordinate in the half
% (upper triangular) Adj matrix
% [coords_i, coords_j] = ind2sub([n, n], find(A));
% [tri_i, tri_j] = ind2sub([n, n], find(triu(A)));

n = size(A, 1);

map = zeros(length(coords_i), 1);
triIdx = containers.Map('KeyType', 'uint32', 'ValueType', 'uint32');

for ind = 1:length(tri_i)
    ii = tri_i(ind);
    jj = tri_j(ind);
    triIdx(n*ii + jj) = ind;
end

for ind = 1:length(coords_i)
    ii = coords_i(ind);
    jj = coords_j(ind);
    if ii < jj
        map(ind) = triIdx(n*ii + jj);
    else
        map(ind) = triIdx(n*jj + ii);
    end
end
