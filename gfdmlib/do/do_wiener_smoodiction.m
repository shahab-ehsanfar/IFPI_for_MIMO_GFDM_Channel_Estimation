function [ H_all_blks, Covariance] = do_wiener_smoodiction(H_p, R_hh_err,pilot_blk, no_p_blk, noblk, N, Ncp, n_H, Fs, fd, numPaths)
% Wiener filtering for joint Smoothing & Prediction for time-variant
% channels over multiple blocks
%
% H_p: Estimated channel frequency repsonse at P-type blocks
% R_hh_err: time domain covariance of channel estimation
% pilot_blk: set of P-type blocks
% no_p_blk: number of P-type blocks
% noblk: total number of blocks
% N: Length of CP-less block (symbol) in samples
% Ncp: Length of CP in samples
% n_H: sample index in which channel estimation is at its best
% Fs: Sampling frequency
% fd: Maximum Doppler shift
% numPaths: Number of channel paths (taps)

H_in_time_with_zeros = ifft(H_p);
H_in_time = H_in_time_with_zeros(1:numPaths,:);

all_blk = 1:noblk;


% n_H: sample index where channel estimation is at its best
N_tot = N + Ncp;
R_f_indx = (pilot_blk-1)*N_tot+ n_H;
R_f_indx_all = (all_blk-1)*N_tot+ n_H;


% Channel time autocorrelation at pilot blocks
R_hp = zeros(no_p_blk);
for i = 1:no_p_blk
    for j = 1:no_p_blk
        R_hp(i,j) = besselj(0,2*pi*( R_f_indx(i) - R_f_indx(j) )/Fs*fd);
    end
end

% Channel time autocorrelation at all blocks
R_h = zeros(noblk);
for i =1:noblk
    for j = 1:noblk
        R_h(i,j) = besselj(0,2*pi*( R_f_indx_all(i) - R_f_indx_all(j) )/Fs*fd);
    end
end

H_smoodicted_time = zeros(numPaths, noblk);
Covariance_l = zeros(noblk,noblk,numPaths);
Covariance = zeros(numPaths,numPaths,noblk);
for l=1:numPaths
   Wiener_mtx = (R_h(pilot_blk,:)'/(R_hp+ R_hh_err(l,l)*eye(no_p_blk)));
   H_smoodicted_time(l,:) = Wiener_mtx*H_in_time(l,:).';
   Covariance_l(:,:,l) = R_h - Wiener_mtx*R_h(pilot_blk,:);   
   for blk = 1:noblk
       Covariance(l,l,blk) = Covariance_l(blk,blk,l);
   end   
end


H_all_blks = fft(H_smoodicted_time,N);




 
end

