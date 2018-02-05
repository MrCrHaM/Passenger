function params = getGraphParams(Adj,p)

n = size(Adj,1);
% unnormalized graph Laplacian
d = sum(Adj,2);
D = diag(d);
L = D-Adj;
% oriented incidence matrix
rows=1;
DI = zeros(1,n);
for i=1:n-1
    for j=i+1:n
        if Adj(i,j)==1
            DI(rows,i)=1;
            DI(rows,j)=-1;
            rows=rows+1;
        end
    end
end
%%%%%%% Adjacency List  %%%%%%%
AL = cell(n,1);
for i=1:n
    AL{i} = find(Adj(i,:)==1);
end

%% compact primal optimization;  fn: node indicator; fe: edge indicator
% pre-processing,  must be done if anchor is changed
M = size(DI,1);         % total #edges
DI_p = DI( abs(DI(:,p))==0 ,:);
params.m = size(DI_p,1);       % #edges after removing node p
DI_star = -eye(n);
DI_star(:,p) = DI_star(:,p)+ones(n,1);  % incidence matrix of p-star graph
Np = zeros(n,1);
Np(find(Adj(p,:)>0)) = 1;
DI_pg = diag(Np)*DI_star;       % incidence matrix of edges with p
params.DI_M = [DI_star; DI_p];             % incidence matrix of M
params.DI_AoM = [DI_pg; DI_p];           % incidence matrix of AoM 
params.DI_p1 = DI_p>0.5;
params.DI_p2 = DI_p<-0.5;
