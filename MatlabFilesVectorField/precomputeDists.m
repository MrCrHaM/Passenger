function d = precomputeDists(tri_i, tri_j, map, y)
% Computes ||y(i,:) - y(j,:)||^2 for any edge (i,j) in graph, returned in
% a m x 1 vector with ordering given by the find(A) function.
% Actually computes above only for upper triangular half of Adj matrix,
% then replicates to lower triangular part using precomputed map.

m = length(tri_i);
d = zeros(m, 1);
for ind = 1:m
    ii = tri_i(ind);
    jj = tri_j(ind);
    d(ind) = sum((y(ii,:) - y(jj,:)).^2);
end

d = d(map);
