close all
clear all

% cvx_solver sedumi

if matlabpool('size') == 0
    matlabpool('open', 12);
end

nSamp = 24;
% gammas = logspace(-2,log10(5),10);
[A, s, yy] = genCounty(2, 1);
pop_vec = ones(1,size(A,1));
S = [69,66,105,93,106,95,49,90,68,59,104,62,88,108,84,92];

params = SAparams(A);
params.l_0 = 0;
% tic;
% [v_Rect, id_Rect] = scan_Rect_RG_real( yy_1,pop_vec, A, params.D_path, params.A_power, params.n, params.l_0);
% 
% [v_SA, id_SA] = SA_fun_real(yy_1, pop_vec,A,params.AL, params.n, 20, S, 0,v_Rect, id_Rect, params.l_0);
% toc;

%%

[A, s, yy] = genCounty(1, nSamp);
scores_noise_rect = zeros(nSamp,1);
scores_noise_SA = zeros(nSamp,1);
parfor ns = 1:nSamp
    ys = yy(:,ns);
    tic;
    [v_Rect, id_Rect] = scan_Rect_RG_real( ys,pop_vec, A, params.D_path, params.A_power, params.n, params.l_0);
    scores_noise_rect(ns) = v_Rect;

    [v_SA, id_SA] = SA_fun_real(ys, pop_vec,A,params.AL, params.n, 40, S, 0,v_Rect, id_Rect, params.l_0);
    scores_noise_SA(ns) = v_SA;
    toc;
end

disp('Noise done')

%%
[A, s, yy] = genCounty(1.5, nSamp);
scores_signal5_rect = zeros(nSamp,1);
scores_signal5_SA = zeros(nSamp,1);
parfor ns = 1:nSamp
    ys = yy(:,ns);
    tic;
    [v_Rect, id_Rect] = scan_Rect_RG_real( ys,pop_vec, A, params.D_path, params.A_power, params.n, params.l_0);
    scores_signal5_rect(ns) = v_Rect;

    [v_SA, id_SA] = SA_fun_real(ys, pop_vec,A,params.AL, params.n, 40, S, 0,v_Rect, id_Rect, params.l_0);
    scores_signal5_SA(ns) = v_SA;
    toc;
end

disp('Signal 1.5 done')

%%
[A, s, yy] = genCounty(1.3, nSamp);
scores_signal3_rect = zeros(nSamp,1);
scores_signal3_SA = zeros(nSamp,1);
parfor ns = 1:nSamp
    ys = yy(:,ns);
    tic;
    [v_Rect, id_Rect] = scan_Rect_RG_real( ys,pop_vec, A, params.D_path, params.A_power, params.n, params.l_0);
    scores_signal3_rect(ns) = v_Rect;

    [v_SA, id_SA] = SA_fun_real(ys, pop_vec,A,params.AL, params.n, 40, S, 0,v_Rect, id_Rect, params.l_0);
    scores_signal3_SA(ns) = v_SA;
    toc;
end

disp('Signal 1.3 done')

%%
[A, s, yy] = genCounty(1.1, nSamp);
scores_signal1_rect = zeros(nSamp,1);
scores_signal1_SA = zeros(nSamp,1);
parfor ns = 1:nSamp
    ys = yy(:,ns);
    tic;
    [v_Rect, id_Rect] = scan_Rect_RG_real( ys,pop_vec, A, params.D_path, params.A_power, params.n, params.l_0);
    scores_signal1_rect(ns) = v_Rect;

    [v_SA, id_SA] = SA_fun_real(ys, pop_vec,A,params.AL, params.n, 40, S, 0,v_Rect, id_Rect, params.l_0);
    scores_signal1_SA(ns) = v_SA;
    toc;
end

disp('Signal 1.1 done')

%%

save county_auc_SA_nopop_36
%%

aucs = computeAUC(scores_noise_rect, scores_signal5_rect, 1);
figure, semilogx(gammas, aucs,'x-')

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