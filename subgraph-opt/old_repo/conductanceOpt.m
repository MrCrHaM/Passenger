function [M, resistance] = conductanceOpt(A, C, gamma, s)

n = size(A,1);

max_iter = 500;
eta = 10;

tic
X = zeros(n,n,max_iter);
temp = ones(n);
X(:,:,1) = temp/trace(temp);
% Xm = X(:,:,1);
Gm = zeros(n,n)/n;

for i = 2:max_iter
    A_M = A.*X(:,:,i-1);
    L_X = diag(sum(A_M,2)) - A_M;
    m = -diag(X(:,:,i-1));
    m(s) = -sum(m([1:s-1,s+1:end]));
    %v = L_X\m %v = v - v(s);
%     v = bicgstab(L_X, m);
    v = pinv(L_X)*m;


    V = repmat(v, 1, n);
    temp = zeros(n); temp(s,s) = v(s)^2/gamma;
    Gc = (2*diag(-v(s) + v) + ((V - V').^2).*A + temp);
    e1 = eigs(Gc,1,'la');
    if m'*v < 1/gamma
        y = 0;%p/20;
    else
        y = norm(C,'fro') / e1;
    end
%     disp(y)
    gradient = C + y * Gc;
    Gm = Gm + gradient; %Gm*(i-1)/i + gradient/i;

    L = eigs(Gm,1,'la');
    X(:,:,i) = expm(eta*(Gm - L*eye(n)));
    X(:,:,i) = X(:,:,i)/trace(X(:,:,i));

%     Xm = Xm*(i-1)/i + X(:,:,i)/i;
end
toc

M = X(:,:,end);
resistance = m'*v;