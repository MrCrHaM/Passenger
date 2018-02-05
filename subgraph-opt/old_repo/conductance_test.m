clear all;
close all;
%cvx_solver sedumi

%% data generation, generate n1xn2 lattice, 4-connected
% n1 = 5; n2 = 5;  % smaller example
%n1 = 5; n2 = 8;  % larger example
n1 = 5; n2 = 1;  % line example
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

% smaller example
% xx_line = zeros(n,1);
% xx_line([7, 9, 10, 15]) = 1;  % 4 nodes with value 1, 1 node disconnected
% xx = reshape(xx_line, n1, n2);
% s = 9;

% larger example
xx(2:3, 2:4) = 1;  % larger connected component
xx(4, 6:8) = 1;  % smaller connected component
s = 13;

% line example
% xx([2:3 6:8]) = 1;
% K = 6;

%rng(1)
% observations with some noise
yy = xx;% + poissrnd(10, size(xx))/100;
% yy(yy < 0) = 0;

% plot using imagesc to display as matrix
% figure, imagesc(yy), title('Original'), colormap gray, colorbar, axis image;

%% enumerate all possible patterns
patts = dec2bin(0:2^n-1)-'0';
%patts(2:end,:) = patts(2:end,:)./repmat(sum(patts(2:end,:),2),1,n);

%warning off
for p = 2:size(patts,1)
    pattern = patts(p,:);
    disp('')
    disp('Pattern:'), disp(pattern)
    M = pattern'*pattern/sum(pattern);
    A_M = Adj.*M;
    A_sub = Adj(find(pattern),find(pattern));
    L_M = diag(sum(A_M,2)) - A_M;
    connected = sum(pattern) == 1 || eigs(diag(sum(A_sub))-A_sub,2,'sm')(2) > eps;
    for s = 1:n
        m = -diag(M);
        m(s) = -sum(m([1:s-1,s+1:end]));
        v = L_M\m;%bicgstab(L_M,m);
        printf('  root: %d\t resistance: %.2f\n', s, m'*v)
    end
end
