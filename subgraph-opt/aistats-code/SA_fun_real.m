function [v_SA, id_SA] = SA_fun_real(yy,pop_vec, Adj,AL, n, N, S,Size_constraint,v_Rect, id_Rect, l_0)
% simulated annealing

% Yuting Chen
% Modified May 26th 2013

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% find the approximate maximum likelihood cluster %%%%
%  yy:          noisy data, column vector, of size nx1
%  Adj:         unweighted adjacency matrix (undirected)
%  n:           size of data
%  k:           searching size parameter, set to be the size of ground
%               truth cluster
%  N:           # of restart the routine
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  v:             optimal cost of indicators: sum(xx(id))/sqrt( sum(id) )
%  id:            logical indicator column vector:  nx1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% n = 24;
% SNR = 4;
% Adj = [0,1,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;1,0,1,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;0,1,0,1,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;0,0,1,0,1,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0;0,0,0,1,0,1,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0;0,0,0,0,1,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0;1,0,0,0,0,0,0,1,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0;0,1,0,0,0,0,1,0,1,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0;0,0,1,0,0,0,0,1,0,1,0,0,0,0,1,0,0,0,0,0,0,0,0,0;0,0,0,1,0,0,0,0,1,0,1,0,0,0,0,1,0,0,0,0,0,0,0,0;0,0,0,0,1,0,0,0,0,1,0,1,0,0,0,0,1,0,0,0,0,0,0,0;0,0,0,0,0,1,0,0,0,0,1,0,0,0,0,0,0,1,0,0,0,0,0,0;0,0,0,0,0,0,1,0,0,0,0,0,0,1,0,0,0,0,1,0,0,0,0,0;0,0,0,0,0,0,0,1,0,0,0,0,1,0,1,0,0,0,0,1,0,0,0,0;0,0,0,0,0,0,0,0,1,0,0,0,0,1,0,1,0,0,0,0,1,0,0,0;0,0,0,0,0,0,0,0,0,1,0,0,0,0,1,0,1,0,0,0,0,1,0,0;0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,1,0,1,0,0,0,0,1,0;0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,1,0,0,0,0,0,0,1;0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,1,0,0,0,0;0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,1,0,1,0,0,0;0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,1,0,1,0,0;0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,1,0,1,0;0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,1,0,1;0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,1,0];
% S = [1,2,3,9,10,16,17,23,24];
% mu = SNR/sqrt(length(S));
% xx_line = zeros(n,1);
% xx_line(S) = mu;
% n1 = 4;
% n2 = 6;
% k = length(S);
% yy = randn(n, 1)+xx_line;
max_iter = 150;
% N = 20; % restart the process N times
vb_th = 10;
%Size_constraint = length(S);

%%%%%%% Adjacency List  %%%%%%%
% AL = cell(n,1);
% for i=1:n
%     AL{i} = find(Adj(i,:)==1);
% end

% likelihood function definition
%likelihood = @(x) sum(yy(x))/sum(pop_vec(x));
 likelihood = @(x) ( (sum(yy(x)) / sum(pop_vec(x)) - l_0 )* sqrt(sum(pop_vec(x))) ) ;

% use the output of rectangle scan as starting point
%[ v_Rect  id_Rect] = scan_Rec( yy, n1, n2 ); 
% IP: just to compare with SA
%[v_IP id_IP IP_gap] = IP_EC( yy, Adj, n, C1, K, ep_IP, S_size,5000,0.08, 0);

% This records the sequence of cluster candidates
S_seq = cell(N, max_iter);
% hL(t, i) = 1 means there exists a neighbor with higher L-value at the current step
hL = logical(zeros(N, max_iter));
% cs(t, i) is the number of consecutive steps such that no new candidates
% with L-value > 1
cs = sparse(zeros(N, max_iter));
% vb(t, i) is the number of times that the current subgraph has been
% visited before in the t-th run
vb = sparse(zeros(N, max_iter));
% cv(t, i) is the common vertices between the current subgraph and the
% best subgraph yet
cv = sparse(zeros(N, max_iter));
% v(t, i) is the value of the likelihood function for i-th iteration at
% t-th run
v = zeros(N, max_iter);
Temp = 'L'; % starting with a high temperature
NewAdd = []; % this is the node that is recently added

shrinking = false;
%tic

for t = 1:N
    shrinking = false;
    if t == 1
        S_seq{t, 1} = find(id_Rect)';
        v(t, 1) = v_Rect;
    else
        pre = [S_seq{2:t-1, 1}]; % the nodes that already used as starting nodes
        % we want to choose a different starting point each time
        this_node =  randi([1, n]);
        while ismember(this_node, pre)
            this_node = randi([1, n]);
        end
        S_seq{t, 1} = this_node;
        v(t,1) = likelihood(this_node);
        NewAdd = [];
    end

    % At each iteration add to remove a node 
    for i = 2:max_iter
        if vb(t, i-1) < vb_th;
        % get all neighboring nodes for the current cluster
        %%%%%%%%%% what about no further operation will increase your value
        to_add = unique(setdiff([AL{[S_seq{t, i-1}]}],[S_seq{t, i-1}]));
        
        %%% begin shrinking since size is bigger than required
        if Size_constraint & vb(t,i-1) >=3 & length(S_seq{t,i-1}) > Size_constraint
            shrinking = true;
        end
        
        if shrinking & length(S_seq{t,i-1}) <= Size_constraint
            shrinking = false;
        end
        if shrinking
            to_add = [];          
        end
        
        
        % get all nodes that could be removed for the current cluster
        % you want to check connectivity for each choice
        if length(S_seq{t, i-1}) > 1
            to_remove = [S_seq{t, i-1}];
        else
            to_remove = [];
        end
        ttt = zeros(1, length(to_remove));
        for jj = 1:length(to_remove)
            S_temp = setdiff(S_seq{t, i-1}, to_remove(jj));
            if length(unique(components(sparse(Adj(S_temp,S_temp))))) > 1
                ttt(jj) = 1; % meaning we'll exclude it from the list
            end
        end
        to_remove(ttt>0) = [];
        %%%%%%%%%% produce different candidate according to temperature
        
        if length(to_add)+length(to_remove) > 0
            
        
        % candidates of clusters at this iteration
        S_list = cell(1, length(to_add)+length(to_remove));
        for jj = 1:length(to_add)
            S_list{1,jj} = unique([S_seq{t, i-1} to_add(jj)]);
        end
        for jj = 1:length(to_remove)
            S_list{1, jj+length(to_add)} = setdiff(S_seq{t, i-1}, to_remove(jj));
        end
        v_list = cellfun(likelihood, S_list);
        [m_v in_v] = max(v_list);
        
        can_list = [to_add to_remove];
        
        else
            break;
            
        end
        
        
        switch Temp
            case 'H' % high temperature, complete random
                p = randperm(length(to_add)+length(to_remove));
                
                in_v = (p(1));
                v(t, i) = v_list(in_v);
                S_seq{t, i} = S_list{in_v};
                for jj = 1:i-1
                    if length(S_seq{t, i}) == length(S_seq{t, jj}) & S_seq{t, i} == S_seq{t, jj}
                        vb(t, i) = vb(t, i) + 1;
                    end
                end
                
                if in_v <= length(to_add)
                    NewAdd = can_list(in_v);
                else
                    NewAdd = [];
                end
                
                if v(t,i) > v(t, i-1)
                    hL(t,i-1) = 1;
                    cs(t, i) = 0;
                else
                    cs(t, i) = cs(t, i-1) + 1;
                end
                
            case 'M' % medium temperature, random in proportion to the likelihood value
                %v_list_norm = v_list./{}
                good_in = find(v_list>0);
                v_list_norm = v_list(good_in)./sum(v_list(good_in));
                v_sum = cumsum(v_list_norm);
                choice = rand;
                for jj = 1:length(good_in)
                    if choice <= v_sum(jj)
                        in_v = (good_in(jj));
                        S_seq{t,i} = S_list{good_in(jj)};
                        v(t, i) = v_list(in_v);
                        break;
                    end
                
                end
                
                if in_v <= length(to_add)
                    NewAdd = can_list(in_v);
                else
                    NewAdd = [];
                end
                
                if v(t,i) > v(t, i-1)
                    hL(t,i-1) = 1;
                    cs(t, i) = 0;
                else
                    cs(t, i) = cs(t, i-1) + 1;
                end
                
                for jj = 1:i-1
                    if length(S_seq{t, i}) == length(S_seq{t, jj}) & S_seq{t, i} == S_seq{t, jj}
                        vb(t, i) = vb(t, i) + 1;
                    end
                end
                
            case 'L' % low temperature, choose the best neighbor
                
               
                v(t, i) = m_v;
                S_seq{t, i} = S_list{in_v};
                
                for jj = 1:i-1
                    if length(S_seq{t, i}) == length(S_seq{t, jj}) & S_seq{t, i} == S_seq{t, jj}
                        vb(t, i) = vb(t, i) + 1;
                    end
                end
                if in_v <= length(to_add)
                    NewAdd = can_list(in_v);
                else
                    NewAdd = [];
                end
                if m_v > v(t, i-1)
                    hL(t,i-1) = 1;
                    cs(t, i) = 0;
                else
                    cs(t, i) = cs(t, i-1) + 1;
                end
                
            case 'C' % call process H(G), continute further in the current direction
                % if a vertice recently added, add one of its neighbours
                % if a vertice recently removed, choose the best neighbor
                if NewAdd
                % means a node is added in last iteration
                    to_add_v = setdiff(AL{NewAdd},S_seq{t,i-1});
                    % this means we do have some neighbor of NewAdd to add
                    if length(to_add_v) > 0
                        p = randperm(length(to_add_v));
                        NewAdd = to_add_v(p(1));
                        S_seq{t, i} = unique([S_seq{t, i-1} NewAdd]);
                        v(t, i) = likelihood(S_seq{t, i});
                        if v(t,i) > v(t, i-1)
                            hL(t,i-1) = 1;
                            cs(t, i) = 0;

                        else
                            cs(t, i) = cs(t, i-1) + 1;
                        end
                        for jj = 1:i-1
                        if length(S_seq{t, i}) == length(S_seq{t, jj}) & S_seq{t, i} == S_seq{t, jj}
                            vb(t, i) = vb(t, i) + 1;
                        end
                        end
                        % otherwise we'll just use temp - L
                    else
                        
                        v(t, i) = m_v;
                        S_seq{t, i} = S_list{in_v};
                        for jj = 1:i-1
                            if length(S_seq{t, i}) == length(S_seq{t, jj}) & S_seq{t, i} == S_seq{t, jj}
                                vb(t, i) = vb(t, i) + 1;
                            end
                        end
                        if in_v <= length(to_add)
                            NewAdd = can_list(in_v);
                        else
                            NewAdd = [];
                        end
                        if m_v > v(t, i-1)
                            hL(t,i-1) = 1;
                            cs(t, i) = 0;
                        else
                            cs(t, i) = cs(t, i-1) + 1;
                        end
                    end 
                 
                 
                else
                    v(t, i) = m_v;
                    S_seq{t, i} = S_list{in_v};
                    for jj = 1:i-1
                        if length(S_seq{t, i}) == length(S_seq{t, jj}) & S_seq{t, i} == S_seq{t, jj}
                            vb(t, i) = vb(t, i) + 1;
                        end
                    end
                    if in_v <= length(to_add)
                        NewAdd = can_list(in_v);
                    else
                        NewAdd = [];
                    end
                    if m_v > v(t, i-1)
                        hL(t,i-1) = 1;
                        cs(t, i) = 0;
                    else
                        cs(t, i) = cs(t, i-1) + 1;
                    end
                    
                end
                
        end
            
        %%%%%%%%%% decide the temperature for next iteration
        %  v(t,i)-likelihood(S) < -1 means the likelihood value is
        %  relatively low
        like_th = 0.2;
        if v(t,i)-likelihood(S) < -like_th  & vb(t,i) >= 5 & cs(t, i) >= 4
            Temp = 'H';
        elseif v(t,i)-likelihood(S) < -like_th & vb(t,i) >=5 & sum(hL(t,max(i-4,1):i))>0
            Temp = 'M';
        elseif  sum(hL(t,max(i-4,1):i))>0 & (v(t,i)-likelihood(S) < -like_th |  vb(t,i) >= 5)
            Temp = 'L';
        elseif v(t,i)-likelihood(S) > -like_th & vb(t,i) < 5 & sum(hL(t,max(i-4,1):i))>0
            Temp = 'C';         
        end
       
        %Temp



        else
            break;
        end
    end

end

% the subgraph that produces the largest value is chosen as the output
% cluster
v_values = cellfun(likelihood, S_seq);
if ~Size_constraint
[v_SA in_v] = max(v_values(:));
f_SA = S_seq{in_v};
id_SA = zeros(n, 1);
id_SA(f_SA) = 1;
else
    le = cellfun(@length, S_seq);
    le_con = find(le<=Size_constraint );
    [v_SA in_v] = max(v_values(le_con));
    f_SA = S_seq{le_con(in_v)};
    id_SA = zeros(n, 1);
    id_SA(f_SA) = 1;
    
end

%v_SA
%v_Rect

%toc



% 
% function re = likelihood(S_can)
% %likelihood = @(x) sum(yy(x))/sum(pop_vec(x));
% 
% if sum(yy(S_can))/sum(pop_vec(S_can)) > sum(yy(setdiff(1:n,S_can)))/sum(pop_vec(setdiff(1:n,S_can)))
% 
% re = sum(yy(S_can)) * log(sum(yy(S_can)) / sum(pop_vec(S_can))) ...
%     + (sum(pop_vec(S_can))-sum(yy(S_can))) * log( (sum(pop_vec(S_can))- sum(yy(S_can))) / sum(pop_vec(S_can))) ...
%     + sum(yy(setdiff(1:n, S_can))) * log( (sum(yy(setdiff(1:n, S_can)))) / (sum(pop_vec(setdiff(1:n, S_can)))) ) ...
%     + ( sum(pop_vec(setdiff(1:n, S_can))) - sum(yy(setdiff(1:n, S_can))) ) * log( ( sum(pop_vec(setdiff(1:n, S_can))) - sum(yy(setdiff(1:n, S_can))) )  / sum(pop_vec(setdiff(1:n, S_can))) );
% 
% else
%     re = sum(yy) * log(sum(yy)) + (sum(pop_vec)-sum(yy)) * log(sum(pop_vec)-sum(yy)) ...
%         - sum(pop_vec) * log(sum(pop_vec));
% 
% end
% 
% % re = ( sum(yy(S_can)) / sum(pop_vec(S_can)) ) * sqrt(sum(pop_vec(S_can)))
% 
% end

end

