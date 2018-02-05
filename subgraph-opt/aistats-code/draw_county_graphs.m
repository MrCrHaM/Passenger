xval = 4;
nSamp = 1;
load('county/map_two_clusters_final.mat')
S = [69,66,105,93,106,95,49,90,68,59,104,62,88,108,84,92];
n = length(cons);

l_0 = 5e-5;
lam_0 = pop_vec * l_0;
l_1 = xval * l_0;
lam_1 = lam_0;
lam_1(S) = pop_vec(S) * l_1;
yy_0 = poissrnd(repmat(lam_0,1,nSamp)); % null hypothesis
yy_1 = poissrnd(repmat(lam_1,1,nSamp));

Adj = zeros(n);
Lats = {cons.Lat};
Lons = {cons.Lon};
fLat = zeros(n,1);
fLon = zeros(n,1);
for i = 1:n
    cons(i).rate = yy_1(i,:)/cons(i).pop;
    cons(i).diseased = 0.5*double(ismember(i, S));
    cons(i).idx = i/n;
    fLat(i) = nanmean(cons(i).Lat);
    fLon(i) = nanmean(cons(i).Lon);
    for j = i+1:n
        Adj(i,j) = double(any(ismembertol(Lats{i},Lats{j},1e-10) & ismembertol(Lons{i},Lons{j},1e-10)));
    end
end

Adj(9,12) = 1;
Adj(12,18) = 1;

Adj = Adj + Adj';
% figure, spy(Adj);
% figure, gplot(Adj, [fLon, fLat],'-*'), text(fLon, fLat+0.1, cellstr(num2str([1:n]')));
s = 90;
yy = reshape([cons.rate],nSamp,n)';
% yy = yy/max(yy);
cons(s).diseased = 1;

fall = flipud(autumn(numel(cons)));
minpop = min([cons.rate]);
maxpop = max([cons.rate]);
densityColors = makesymbolspec('Polygon', {'rate', [0 maxpop], 'FaceColor', fall});
figure, geoshow(cons, 'DisplayType', 'polygon','SymbolSpec', densityColors),
caxis([0 maxpop]), colormap(fall), set(gca,'pos',[0 0 0.9 0.9]), axis off, colorbar

cm = [fall(1,:); fall(round(end/2),:); fall(end,:)];
diseasedColors = makesymbolspec('Polygon', {'diseased', [0 1], 'FaceColor', cm});
figure, geoshow(cons, 'DisplayType', 'polygon','SymbolSpec', diseasedColors), colormap(cm);
set(gca,'pos',[0 0 0.9 0.9]), axis off;
h = lcolorbar({'Normal', 'Anomalous', 'Anchor'}); set(h,'position',[0.8 0.2 0.04 0.6]);

G = graph(Adj);
colors = zeros(n,1);
colors(S) = 0.5;
colors(s) = 1;
figure, plot(G, 'XData', fLon, 'YData', fLat, 'NodeCData', colors, 'MarkerSize', 4), colormap(fall), axis image, axis off