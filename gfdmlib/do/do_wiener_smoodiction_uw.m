function [ h_all_blks, Covariance] = do_wiener_smoodiction_uw(h_hat, R_hh_err, noblk, N, n_H, n_Hd, Fs, fd, numPaths)
% Wiener filtering for joint Smoothing & Prediction for time-variant
% channels over multiple blocks
%
% Shahab Ehsanfar, TU Dresden
% n_H: sample index where channel estimation is at its best

H_in_time = h_hat;
no_p_blk = size(H_in_time,2);

all_blk = 1:noblk;
pilot_blk = 1:no_p_blk;

% R_hh_err = zeros(size(R_hh_err));
N_tot = N;
R_f_indx = (pilot_blk-1)*N_tot+ n_H;
R_f_indx_all = (all_blk-1)*N_tot+ n_H;

R_f_indx_all = [R_f_indx_all(1:end-1) (noblk-2)*N_tot+n_Hd R_f_indx_all(end)];
pilot_blk = [1:no_p_blk-1 no_p_blk+1];

% Channel time autocorrelation at pilot blocks
R_hp = zeros(no_p_blk);
for i = 1:no_p_blk
    for j = 1:no_p_blk
        R_hp(i,j) = besselj(0,2*pi*( R_f_indx(i) - R_f_indx(j) )/Fs*fd);
    end
end

% Channel time autocorrelation at all blocks
no_wblk = noblk+1;
R_h = zeros(no_wblk);
for i =1:no_wblk
    for j = 1:no_wblk
        R_h(i,j) = besselj(0,2*pi*( R_f_indx_all(i) - R_f_indx_all(j) )/Fs*fd);
    end
end

H_smoodicted_time = zeros(numPaths, no_wblk);
Covariance_l = zeros(no_wblk,no_wblk,numPaths);
Covariance = zeros(numPaths,numPaths,no_wblk);
for l=1:numPaths
   Wiener_mtx = (R_h(pilot_blk,:)'/(R_hp+ R_hh_err(l,l)*eye(no_p_blk)));
   H_smoodicted_time(l,:) = Wiener_mtx*H_in_time(l,:).';
   Covariance_l(:,:,l) = R_h - Wiener_mtx*R_h(pilot_blk,:);
   % mse(l) = abs(mean(diag(Covariance(:,:,l))));
   for blk = 1:no_wblk
       Covariance(l,l,blk) = Covariance_l(blk,blk,l);
   end   
end
% figure; plot(mse)

h_all_blks = H_smoodicted_time; 

if noblk == 19
   []; 
end

%Covariance = mean(Covariance_l,3);
 
end

