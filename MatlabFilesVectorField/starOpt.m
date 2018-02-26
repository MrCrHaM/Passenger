function [M,Y] = starOpt(A, C, gamma, s)

n = size(A,1);
degs = sum(A)'; 
% degs = ones(n, 1);

% [~,s] = max(diag(C));

max_iter = 100;
eta = 5;
p = 2*norm(C,'fro')/gamma;

tic
% randomly init Y
Y = p/sum(degs)*eye(n); %(eye(n)+0.1*rand(n)); Y = Y'*Y; Y = Y/trace(Y)*p;
X = zeros(n,n,max_iter);
grad = zeros(n);

opts.issym = 1;
for i = 1:max_iter
    % first eigenvector matrix for primal opt
    [u,~] = eigs(real(paramQ_star(A, degs, gamma, Y, s) + C), 1, 'la', opts);
    opts.v0 = u;
    X(:,:,i) = u*u';
    QX = paramQ_star(A, degs, gamma, [], s, X(:,:,i));
    
    % MWU for dual variable Y
    grad = grad - QX;
    Y = expm(eta*grad/max(degs)); %logm(Y) - eta*QX);
    
    Y = Y/trace(diag(degs)*Y)*p;  % normalize to trace p
end
M = mean(X,3);
toc

% Qf = paramQ_11(Adj,gamma,[],r,final);
% ee = min(eig(Qf));

% figure, imagesc(reshape(diag(final), n1, n2)), colorbar, title('Estimated');
% figure, imagesc(final), colorbar

