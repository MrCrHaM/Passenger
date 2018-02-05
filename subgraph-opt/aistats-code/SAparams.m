function params = SAparams(Adj)

n = size(Adj,1);
params.D_path = all_shortest_paths(sparse(full(Adj)));
% max_dist = max(D_path(:));
max_dist =  max(params.D_path(params.D_path<Inf));

% A_power = cell(1, max_dist);
% for i = 1:max_dist
%     A_power{i} = full(logical(Adj^i));
% end
params.A_power = (zeros(n,n,max_dist));
for i = 1:max_dist
    params.A_power(:,:,i) = full(logical(Adj^i));
end
params.A_power = logical(params.A_power);

%%%%%%% Adjacency List  %%%%%%%
params.AL = cell(n,1);
for i=1:n
    params.AL{i} = find(Adj(i,:)==1);
end

% K = max(cellfun('length', AL))+2;          % max degree parameter (is actually Max degree + 2)

params.l_0 = 5e-5;
params.n = n;