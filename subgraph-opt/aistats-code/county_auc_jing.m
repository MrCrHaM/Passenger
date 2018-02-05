close all
clear all

cvx_solver sedumi

if matlabpool('size') == 0
    matlabpool('open', 4);
end

nSamp = 20;
gammas = logspace(-2,log10(5),10);
% [A, s, yy] = genCounty(1.1, nSamp);

%%

[A, s, yy] = genCounty(1, nSamp);
params = getGraphParams(A,s);
scores_noise = zeros(nSamp,length(gammas));
for ns = 1:nSamp
    ys = yy(:,ns);
    c = ys/norm(ys)*sqrt(length(ys));
%     C = c*c';
    for gind = 1:length(gammas)
        fn = jing_SDP(A, params, c, 16, gammas(gind), s);
        scores_noise(ns,gind) = ys'*fn;
    end
end

disp('Noise done')

%%
[A, s, yy] = genCounty(1.5, nSamp);
scores_signal5 = zeros(nSamp,length(gammas));
for ns = 1:nSamp
    ys = yy(:,ns);
    c = ys/norm(ys)*sqrt(length(ys));
%     C = c*c';
    for gind = 1:length(gammas)
        fn = jing_SDP(A, params, c, 16, gammas(gind), s);
        scores_signal5(ns,gind) = ys'*fn;
    end
end

disp('Signal 1.5 done')

%%
[A, s, yy] = genCounty(1.3, nSamp);
scores_signal3 = zeros(nSamp,length(gammas));
for ns = 1:nSamp
    ys = yy(:,ns);
    c = ys/norm(ys)*sqrt(length(ys));
%     C = c*c';
    for gind = 1:length(gammas)
        fn = jing_SDP(A, params, c, 16, gammas(gind), s);
        scores_signal3(ns,gind) = ys'*fn;
    end
end

disp('Signal 1.3 done')

%%
[A, s, yy] = genCounty(1.1, nSamp);
scores_signal1 = zeros(nSamp,length(gammas));
for ns = 1:nSamp
    ys = yy(:,ns);
    c = ys/norm(ys)*sqrt(length(ys));
%     C = c*c';
    for gind = 1:length(gammas)
        fn = jing_SDP(A, params, c, 16, gammas(gind), s);
        scores_signal1(ns,gind) = ys'*fn;
    end
end

disp('Signal 1.1 done')

%%

save county_auc_jing
%%

% aucs = computeAUC(scores_noise, scores_signal3, gammas);
% figure, semilogx(gammas, aucs,'x-')

aucs1 = computeAUC(scores_noise, scores_signal1, gammas);
aucs3 = computeAUC(scores_noise, scores_signal3, gammas);
aucs5 = computeAUC(scores_noise, scores_signal5, gammas);
figure, semilogx(gammas, aucs1, 'x-', gammas, aucs3, 'o--', gammas, aucs5, '+:')
legend('\lambda_1/\lambda_0 = 1.1', '\lambda_1/\lambda_0 = 1.3', '\lambda_1/\lambda_0 = 1.5')
xlabel('\gamma'), ylabel('AUC')

% figure, plot(fps,tps,'x-')
%%
% gamma = 0.1;
% 
% ys = yy(:,10);
% c = ys/norm(ys)*sqrt(length(ys));
% C = c*c';
% 
% M = starOpt(A, C, gamma, s);
% 
% figure, 
% subplot(3,1,1), stem(diag(M))
% % subplot(2,1,1), imagesc(reshape(M(s,:), n1, n2)), colormap gray, colorbar, axis image
% subplot(3,1,2), stem(ys)
% subplot(3,1,3), stem(ismember(1:129,[69,66,105,93,106,95,49,90,68,59,104,62,88,108,84,92])); %,80,85,76,107,54,86,71]))
% % subplot(2,1,2), imagesc(reshape(diag(M), n1, n2)), colormap gray, colorbar, axis image,
% % subplot(3,1,3), imagesc(reshape(diag(Gm), n1, n2)), colormap gray, colorbar, axis image,