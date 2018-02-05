function [Adj, pts] = readGeoGraph(filename)

params = dlmread(filename, '', [0 0 0 2]);
n = params(1); d = params(2); m = params(3);
pts = dlmread(filename, '', [1 0 n d-1]);
edges = dlmread(filename, '', [n+1 0 n+m 1]);

Adj = sparse(edges(:,1) + 1, edges(:,2) + 1, ones(m,1), n, n);
Adj = Adj + Adj';