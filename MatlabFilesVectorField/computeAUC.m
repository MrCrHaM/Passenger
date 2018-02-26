function aucs = computeAUC(scores_noise, scores_signal, gammas)

gamma_pts = 500;
aucs = zeros(length(gammas),1);
% figure,
for gind = 1:length(gammas)
    th_range = linspace(mean(scores_noise(:,gind))/10, mean(scores_signal(:,gind))*10, gamma_pts);
    tps = zeros(gamma_pts,1); fps = zeros(gamma_pts,1);
    for thind = 1:gamma_pts
        th = th_range(thind);
        tps(thind) = mean(scores_signal(:,gind) > th);
        fps(thind) = mean(scores_noise(:,gind) > th);
    end
    fps = [fps' 0 1 1]';
    tps = [tps' 0 1 0]';

    K = convhull([fps; 1],[tps; 0]);
    fps = fps(K);
    tps = tps(K);

    auc = polyarea(fps,tps); %trapz(fps_s,tps_s);
    aucs(gind) = auc;
%     plot(fps,tps,'x-'), title(num2str(gammas(gind)))
%     waitforbuttonpress;
end
% figure, semilogx(gammas, aucs,'x-')