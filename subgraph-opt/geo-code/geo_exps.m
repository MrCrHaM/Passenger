close all
clear all

p = gcp('nocreate'); % If no pool, do not create new one.
if isempty(p)
    parpool('local', 4);
end

% generate graph with "geogen 3 500 5 1 > example.graph"

[A, pts] = readGeoGraph('../example_10k_3d.graph');

gammas = logspace(-3,log10(2),10);
nSamp = 50;

%%
[yy, S] = genMeasurements(pts, 40, 100, 100, nSamp);
s = S(1);
scores_noise = zeros(nSamp,10);
parfor ns = 1:nSamp
    ys = yy(:,ns);
    c = ys/norm(ys)*sqrt(length(ys));
%     C = c*c';
    for gind = 1:10
        u_s = starOpt_fast(A, c, gammas(gind), s, 10, 100);
        scores_noise(ns,gind) = mean((ys'*u_s).^2); %trace(ys*ys'*M);
    end
end

%%
[yy, S] = genMeasurements(pts, 40, 100, 150, nSamp);
s = S(1);
scores_signal = zeros(nSamp,10);
parfor ns = 1:nSamp
    ys = yy(:,ns);
    c = ys/norm(ys)*sqrt(length(ys));
%     C = c*c';

    for gind = 1:10
        u_s = starOpt_fast(A, c, gammas(gind), s, 10, 100);
        scores_signal(ns,gind) = mean((ys'*u_s).^2);
    end
end

%%
aucs_fast = computeAUC(scores_noise, scores_signal, gammas);
figure, semilogx(gammas, aucs_fast, 'o--')