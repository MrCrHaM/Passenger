close all
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
[A, s, yy] = genCounty(1.3, nSamp);
scores_signal3 = zeros(nSamp,10);
parfor ns = 1:nSamp
    ys = yy(:,ns);
    c = ys/norm(ys)*sqrt(length(ys));
    C = c*c';

    for gind = 1:10
        M = starOpt(A, C, gammas(gind), s);
        scores_signal3(ns,gind) = trace(ys*ys'*M);
    end
end

disp('Signal 1.3 done')

%%
[A, s, yy] = genCounty(1.1, nSamp);
scores_signal1 = zeros(nSamp,10);
parfor ns = 1:nSamp
    ys = yy(:,ns);
    c = ys/norm(ys)*sqrt(length(ys));
    C = c*c';

    for gind = 1:10
        M = starOpt(A, C, gammas(gind), s);
        scores_signal1(ns,gind) = trace(ys*ys'*M);
    end
end

disp('Signal 1.1 done')

%%

save county_auc_gammas_all_50samp_norm_100iter_new
%%

% gamma_pts = 500;
% scores_signal = scores_signal1;
% aucs = zeros(10,1);
% % figure,
% for gind = 1:10
%     th_range = linspace(mean(scores_noise(:,gind))/10, mean(scores_signal(:,gind))*10, gamma_pts);
%     tps = zeros(gamma_pts,1); fps = zeros(gamma_pts,1);
%     for thind = 1:gamma_pts
%         th = th_range(thind);
%         tps(thind) = mean(scores_signal(:,gind) > th);
%         fps(thind) = mean(scores_noise(:,gind) > th);
%     end
%     fps = [fps' 0 1 1]';
%     tps = [tps' 0 1 0]';
% 
%     K = convhull([fps; 1],[tps; 0]);
%     fps = fps(K);
%     tps = tps(K);
% %     [~, si] = sort(fps);
% %     fps_s = fps(si); tps_s = tps(si);
% 
%     auc = polyarea(fps,tps); %trapz(fps_s,tps_s);
%     aucs(gind) = auc;
% %     plot(fps,tps,'x-'), title(num2str(gammas(gind)))
% %     waitforbuttonpress;
% end

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