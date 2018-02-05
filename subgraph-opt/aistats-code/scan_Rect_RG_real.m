function [v id] = scan_Rect_RG_real( yy,pop_vec, Adj, D_path, A_power, n, l_0);  
% rectangle scan for graphs

% possible candidate S set: from each node, include all nodes within k-hop
% distance for every possible k

% we then pick the S with maximum likelihood 

% Yuting Chen
% May 21, 2013

%%%%%%%%%%%%%%%% Input %%%%%%%%%%%%%%%%%%
% yy: noisy input
% Adj: adjacency matrix for graph
% D_path: shortest path distance for each node pair
% A_power: the list of powers of A (logical) (for determining k-hop neighbourhoods)
% n: number of nodes

%%%%%%%%%%%%%%%% Output %%%%%%%%%%%%%%%%%
% v : likelihood value
% id: binary indicator of S
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% record the maximum k-hop distance from each node to all other nodes
D_path(D_path==Inf) =  max(D_path(D_path<Inf));
max_D = max(D_path);

% then we have this many different S candidates
S_candidate = cell(1, sum(max_D));

iter = 1;
for i = 1:n
    for k = 1:max_D(i)
        % find all the nodes that lie within k-hop distance of node i
        S_temp = [i];
        for j = 1:k
            S_temp = union(S_temp, find(A_power(i,:,j)));
        end
        S_candidate{iter} = S_temp;
        iter = iter + 1;
    end
end

% % test connectivity
% t = zeros(1, length(S_candidate));
% for i = 1:length(S_candidate)
%     t(i) = length(unique(components(sparse(Adj(S_candidate{i},S_candidate{i})))));
%      
% end

% find the candidate that has the largest likelihood value
%likelihood = @(x) sum(yy(x))/sum(pop_vec(x));
likelihood = @(x) ( (sum(yy(x)) / sum(pop_vec(x)) - l_0 )* sqrt(sum(pop_vec(x))) ) ;

v_temp = cellfun(likelihood, S_candidate);

[v_m i_m] = max(v_temp);

v = v_m;
id = zeros(n, 1);
id(S_candidate{i_m}) = 1;
