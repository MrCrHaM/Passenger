close all
clear all

load('county/map_two_clusters_final.mat')

n = length(cons);

l_0 = 5e-5;
lam_0 = pop_vec * l_0;
l_1 = 8 * l_0;
lam_1 = lam_0;
lam_1(S) = pop_vec(S) * l_1;
yy_0 = poissrnd(lam_0); % null hypothesis
yy_1 = poissrnd(lam_1);

Adj = zeros(n);
Lats = {cons.Lat};
Lons = {cons.Lon};
names = {cons.NAME};
fLat = zeros(n,1);
fLon = zeros(n,1);
for i = 1:n
    cons(i).rate = yy_1(i)/cons(i).pop;
    cons(i).diseased = double(ismember(i, S));
    cons(i).idx = i/n;
%     fLat(i) = mean(cons(i).Lat, 'omitnan');
%     fLon(i) = mean(cons(i).Lon, 'omitnan');
    [G, ~] = polygeom2(cons(i).Lon(~isnan(cons(i).Lon)), cons(i).Lat(~isnan(cons(i).Lat)));
    fLon(i) = G.x_c; fLat(i) = G.y_c;
    for j = i+1:n
        Adj(i,j) = double(any(ismembertol(Lats{i},Lats{j},1e-1, 'DataScale', 1) ...
            & ismembertol(Lons{i},Lons{j},5e-2, 'DataScale', 1)));
    end
end
Adj(9,18) = 1; Adj(12,18) = 1;  % connect Nantucket to Barnstable and Dukes
Adj = Adj + Adj';
% figure, spy(Adj);
% figure, gplot(Adj, [fLon, fLat],'-*'), axis equal, xlim([-80, -66]), ylim([40 48]);

% figure, geoshow(cons, 'DisplayType', 'polygon','SymbolSpec', makesymbolspec('Polygon', {'idx', [0 1], 'FaceColor', jet(n)}))
% xlim([-80, -66]), ylim([40 48]), text(fLon, fLat, strsplit(num2str(1:n)));


%%
fall = flipud(autumn(numel(cons)));
minpop = min([cons.rate]);
maxpop = max([cons.rate]);
densityColors = makesymbolspec('Polygon', {'rate', [minpop maxpop], 'FaceColor', fall});
diseasedColors = makesymbolspec('Polygon', {'diseased', [0 1], 'FaceColor', fall});
% figure, geoshow(cons, 'DisplayType', 'polygon','SymbolSpec', diseasedColors);
figure, geoshow(cons, 'DisplayType', 'polygon','SymbolSpec', densityColors);
% figure, geoshow(cons, 'DisplayType', 'polygon','SymbolSpec', makesymbolspec('Polygon', {'idx', [0 1], 'FaceColor', jet(n)}));

yy = [cons.rate];
yy = yy/max(yy);

%%
gamma = 2;

c = yy(:);
C = c*c';
s = 93;

[M, res] = conductanceOpt(Adj, C, gamma, s);

for i = 1:n
    cons(i).est = double(M(s,i) > 0.03);
end

figure, 
% subplot(2,1,1), stem(reshape(diag(final), n1, n2))
geoshow(cons, 'DisplayType', 'polygon','SymbolSpec', makesymbolspec('Polygon', {'est', [0 1], 'FaceColor', fall}))

% imagesc(reshape(M(s,:), n1, n2)), colormap gray, colorbar, axis image
% subplot(2,1,2), stem(yy)
% subplot(2,1,2), imagesc(reshape(diag(M), n1, n2)), colormap gray, colorbar, axis image,
% subplot(3,1,3), imagesc(reshape(diag(Gm), n1, n2)), colormap gray, colorbar, axis image,