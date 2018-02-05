function [n1, n2, Adj, s, yy] = genExample(selector, noisy)

switch selector
    case 'line'
        n1 = 10; n2 = 2;  % line example
    case 'small'
        n1 = 5; n2 = 5;  % smaller example
    case 'large'
        n1 = 5; n2 = 8;  % larger example
end

n = n1*n2;

% compute unweighted adjacency matrix, Laplacian matrix
Adj = zeros(n, n);
for i=1:n1
    for j=1:n2
        if j<n2
            Adj(i+(j-1)*n1, i+j*n1) = 1;
            Adj(i+j*n1, i+(j-1)*n1) = 1;
        end
        if i<n1
            Adj(i+(j-1)*n1, i+1+(j-1)*n1) = 1;
            Adj(i+1+(j-1)*n1, i+(j-1)*n1) = 1;
        end
    end
end

xx = zeros(n1, n2);

switch selector
    case 'line'
        xx([2:3 6:8]) = 1;
        s = 3;
    case 'small'
        xx_line = zeros(n,1);
        xx_line([7, 9, 10, 15]) = 1;  % 4 nodes with value 1, 1 node disconnected
        xx = reshape(xx_line, n1, n2);
        s = 9;
    case 'large'
        xx(2:3, 2:4) = 1;  % larger connected component
        xx(4, 6:8) = 1;  % smaller connected component
        s = 13;
end

yy = xx + noisy * poissrnd(10, size(xx))/100;