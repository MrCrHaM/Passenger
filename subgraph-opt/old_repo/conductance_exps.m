close all
clear all

[n1, n2, A, s, yy] = genExample('small', 1);

gamma = 10;

c = yy(:);
C = c*c';

[M, res] = conductanceOpt(A, C, gamma, s);

figure, 
% subplot(2,1,1), stem(reshape(diag(final), n1, n2))
subplot(2,1,1), imagesc(reshape(M(s,:), n1, n2)), colormap gray, colorbar, axis image
% subplot(2,1,2), stem(yy)
subplot(2,1,2), imagesc(reshape(diag(M), n1, n2)), colormap gray, colorbar, axis image,
% subplot(3,1,3), imagesc(reshape(diag(Gm), n1, n2)), colormap gray, colorbar, axis image,