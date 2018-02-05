function ind = ismembertol_custom(A,B,tol)
if size(A,1) == 1
    A = A';
end
if size(B,1) == 1
    B = B';
end
d = pdist2(A,B);
ind = sum(d <= tol*max(abs([A(:);B(:)])), 2) > 0;