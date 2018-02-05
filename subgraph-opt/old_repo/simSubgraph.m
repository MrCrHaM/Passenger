clear all;
close all;
cvx_solver sedumi

% parameter list
gammas = [0.01 0.1 0.25];
sigmas = [0.01 0.1 0.25 0.5];
numGraphs = 1;

%% data generation, generate n1xn2 lattice, 4-connected
% n1 = 5; n2 = 5;  % smaller example
n1 = 8; n2 = 8;  % larger example
n = n1*n2;

% compute unweighted adjacency matrix, Laplacian matrix
Adj = zeros(n, n);
for i=1:n1
    for j=1:n2
        if j<n2
            Adj((i-1)*n2+j, (i-1)*n2+j+1) = 1;
            Adj((i-1)*n2+j+1, (i-1)*n2+j) = 1;
        end
        if i<n1
            Adj((i-1)*n2+j, i*n2+j) = 1;
            Adj(i*n2+j, (i-1)*n2+j) = 1;
        end
    end
end
% unnormalized graph Laplacian
d = sum(Adj,2);
D = diag(d);
L = D-Adj;

% smaller example
% xx_line = zeros(n,1);
% xx_line([7, 9, 10, 15]) = 1;  % 4 nodes with value 1, 1 node disconnected
% xx = reshape(xx_line, n1, n2);

% larger example
xx = zeros(n1,n2);
xx(3:4, 2:4) = 1;  % larger connected component
xx(5, 6:8) = 1;  % smaller connected component
K = 9;

%% generate and solve examples
% input and outputs storage for each graph, sigma, gamma and sparsity constraint existing or not
% inputs = zeros(numGraphs,n1,n2);
observations = zeros(numGraphs, length(sigmas), n1, n2);
outputs = cell(numGraphs, 1);

for gr_ind = 1%1:numGraphs
    % random rectangle
%     ts = randi(n1-2); bs = randi(n1-2);
%     ls = randi(n2-2); rs = randi(n2-2);
%     t = min(ts,bs); b = max(ts,bs);
%     l = min(ls,rs); r = max(ls,rs);
%     xx = zeros(n1, n2);
%     xx(t+1:b+1,l+1:r+1) = 1;
%     K = sum(xx(:));

%     inputs(gr_ind,:,:) = xx;
    obs = zeros(length(sigmas), n1, n2);
    outs1 = zeros(length(sigmas), length(gammas), n, n);
    outs2 = zeros(length(sigmas), length(gammas), n, n);

    for s_ind = 1:length(sigmas)
        sigma = sigmas(s_ind);
        yy = xx + sigma*randn(size(xx));
        obs(s_ind,:,:) = yy;

        for g_ind = 1:length(gammas)
            gamma = gammas(g_ind);

            % plot using imagesc to display as matrix
            % figure, imagesc(yy), title('Observed');

            % cvx primal opt, with sparsity
            fprintf('\n\n\nRunning for graph %d, sigma = %.2f, gamma = %.2f, no sparsity\n', gr_ind, sigma, gamma)
            M1 = subgraphOpt(Adj, yy(:), 0, gamma);

            % cvx primal opt, with sparsity
            fprintf('\n\n\nRunning for graph %d, sigma = %.2f, gamma = %.2f, with sparsity\n', gr_ind, sigma, gamma)
            M2 = subgraphOpt(Adj, yy(:), K, gamma);

            % save results
            outs1(s_ind,g_ind,:,:) = M1;
            outs2(s_ind,g_ind,:,:) = M2;
        end
    end
    observations(gr_ind,:,:,:) = obs;
    outputs{gr_ind} = {outs1, outs2};
end

save results_fixed.mat observations outputs

%%
load results_fixed.mat
outs1 = outputs{1}{1};
outs2 = outputs{1}{2};

% plot one noisy observation
figure, imagesc(squeeze(observations(1,2,:,:))), colorbar, colormap gray, 
h = title('Observation with $\sigma = 0.1$');
set(h,'interpreter','latex')

% display diagonal rounding solutions without K constraint
figure, ind = 0;
for s_ind = 1:length(sigmas)
    sigma = sigmas(s_ind);
    for g_ind = 1:length(gammas)
        ind = ind + 1;
        gamma = gammas(g_ind);
        M1 = squeeze(outs1(s_ind,g_ind,:,:));
%         subgraphRound(M1, n1, n2)
        
        subplot(length(sigmas),length(gammas), ind),
        imagesc(reshape(diag(M1), n1, n2)), colorbar, colormap gray, title(sprintf('\\sigma = %.2f, \\gamma = %.2f', sigma, gamma));
    end
end

% display diagonal rounding solutions with K constraint
figure, ind = 0;
for s_ind = 1:length(sigmas)
    sigma = sigmas(s_ind);
    for g_ind = 1:length(gammas)
        ind = ind + 1;
        gamma = gammas(g_ind);
        M2 = squeeze(outs2(s_ind,g_ind,:,:));
%         subgraphRound(M1, n1, n2)
        
        subplot(length(sigmas),length(gammas), ind),
        imagesc(reshape(diag(M2), n1, n2)), colorbar, colormap gray, axis image
        h = title(sprintf('$\\sigma = %.2f$, $\\gamma = %.2f$', sigma, gamma));
        set(h,'interpreter','latex')
    end
end

% display nonzero submatrix of solution with K constraint
figure, ind = 0;
for s_ind = 1:length(sigmas)
    sigma = sigmas(s_ind);
    for g_ind = 1:length(gammas)
        ind = ind + 1;
        gamma = gammas(g_ind);
        M2 = squeeze(outs2(s_ind,g_ind,:,:));
%         subgraphRound(M1, n1, n2)
        S = diag(M2) > 1e-3;
        Msub = M2(find(S), find(S));  % non-zero submatrix of M
        subplot(length(sigmas),length(gammas), ind),
        imagesc(Msub), colorbar, colormap gray, title(sprintf('\\sigma = %.2f, \\gamma = %.2f', sigma, gamma));
    end
end


[ind1, ind2] = ind2sub([n1 n2], (1:n)');
ind_text = strcat(num2str(ind1), ',', num2str(ind2));

% display hyperplane rounding to two dim.s of solution with K constraint
figure, ind = 0;
for s_ind = 1:length(sigmas)
    sigma = sigmas(s_ind);
    for g_ind = 1:length(gammas)
        ind = ind + 1;
        gamma = gammas(g_ind);
        M = squeeze(outs2(s_ind,g_ind,:,:));
%         subgraphRound(M1, n1, n2)

        % project and display
        S = diag(M) > 1e-3;
        Msub = M(find(S), find(S));  % non-zero submatrix of M
%         Msub = Msub + 1e-3*eye(length(Msub));
%         Msub = M + 0.001*eye(length(Msub));
        
        % lattice coordinates of indices in M
        ind_sub = ind_text(S, :);

        R = randn(size(Msub,1),2);  % random Gaussian vectors
        R = R./repmat(sqrt(sum(R.^2)), size(Msub,1), 1);  % normalize to unit length

        V = Msub*R;
        subplot(length(sigmas),length(gammas), ind), scatter(V(:,1), V(:,2)), hold on, xlim([-0.3 0.3]), ylim([-0.3 0.3]),
        text(V(:,1), V(:,2), ind_sub), title(sprintf('\\sigma = %.2f, \\gamma = %.2f', sigma, gamma))
    end
end

% display hyperplane rounding to two dim.s of solution with K constraint, but project the sqrt of submatrix
figure, ind = 0;
for s_ind = 1:length(sigmas)
    sigma = sigmas(s_ind);
    for g_ind = 1:length(gammas)
        ind = ind + 1;
        gamma = gammas(g_ind);
        M = squeeze(outs2(s_ind,g_ind,:,:));
%         subgraphRound(M1, n1, n2)

        % project and display
        S = diag(M) > 1e-3;
        Msub = M(find(S), find(S));  % non-zero submatrix of M
        [v d] = eig(Msub);
        d(d < 0) = 1e-8;
        Msub = v*sqrt(d)*v';
%         Msub = Msub + 1e-5*eye(length(Msub));
%         Msub = M + 0.001*eye(length(Msub));
        
        % lattice coordinates of indices in M
        ind_sub = ind_text(S, :);

        R = randn(size(Msub,1),2);  % random Gaussian vectors
        R = R./repmat(sqrt(sum(R.^2)), size(Msub,1), 1);  % normalize to unit length

        Mchol = chol(Msub);
        V = Mchol'*R;
        subplot(length(sigmas),length(gammas), ind), scatter(V(:,1), V(:,2)), hold on, xlim([-0.3 0.3]), ylim([-0.3 0.3]),
        text(V(:,1), V(:,2), ind_sub), title(sprintf('\\sigma = %.2f, \\gamma = %.2f', sigma, gamma))
    end
end
