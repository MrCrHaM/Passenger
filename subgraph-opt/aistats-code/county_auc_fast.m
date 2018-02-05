% close all
clear all

p = gcp('nocreate'); % If no pool, do not create new one.
if isempty(p)
    parpool('local', 12);
end

nSamp = 50;
gammas = logspace(-2,log10(5),10);
% [A, s, yy] = genCounty(1.1, nSamp);

%%
[A, s, yy] = genCounty(1, nSamp);
scores_noise = zeros(nSamp,10);
parfor ns = 1:nSamp
    ys = yy(:,ns);
    c = ys/norm(ys)*sqrt(length(ys));
    C = c*c';
    for gind = 1:10
        M = starOpt(A, C, gammas(gind), s);
        scores_noise(ns,gind) = trace(ys*ys'*M);
    end
end

disp('Noise done')

%%
% [A, s, yy] = genCounty(1, nSamp);
scores_noise_fast = zeros(nSamp,10);
parfor ns = 1:nSamp
    ys = yy(:,ns);
    c = ys/norm(ys)*sqrt(length(ys));
    C = c*c';
    for gind = 1:10
        u_s = starOpt_fast(A, c, gammas(gind), s, 20, 100);
        scores_noise_fast(ns,gind) = mean((ys'*u_s).^2); %trace(ys*ys'*M);
    end
end

disp('Noise fast done')

%%
[A, s, yy] = genCounty(1.5, nSamp);
scores_signal5 = zeros(nSamp,10);
parfor ns = 1:nSamp
    ys = yy(:,ns);
    c = ys/norm(ys)*sqrt(length(ys));
    C = c*c';

    for gind = 1:10
        M = starOpt(A, C, gammas(gind), s);
        scores_signal5(ns,gind) = trace(ys*ys'*M);
    end
end

disp('Signal 1.5 done')

%%
%[A, s, yy] = genCounty(1.5, nSamp);
scores_signal5_fast = zeros(nSamp,10);
parfor ns = 1:nSamp
    ys = yy(:,ns);
    c = ys/norm(ys)*sqrt(length(ys));
    C = c*c';

    for gind = 1:10
        u_s = starOpt_fast(A, c, gammas(gind), s, 20, 100);
        scores_signal5_fast(ns,gind) = mean((ys'*u_s).^2);
    end
end

disp('Signal 1.5 fast done')

%%
aucs5 = computeAUC(scores_noise, scores_signal5, gammas);
aucs5_fast = computeAUC(scores_noise_fast, scores_signal5_fast, gammas);
figure, semilogx(gammas, aucs5, 'x-', gammas, aucs5_fast, 'o--')
