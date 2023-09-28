function [ d_hat , R_dd_tilde, equalizer, R_dd_hat] = do_mimo_uw_equalization( config, snr, H_hat, R_hh_err, Y_all)
% Parallel Interference Cancellation for MIMO GFDM Channel Estimation
%
% TU-Dresden
% Shahab Ehsanfar
%

nT = config.nT;
nR = config.nR;
N = config.N;
Np = config.preamble.N_preamble/2;
N_data = config.payload.N_data;
numPaths = config.numPaths;
FA = config.matrices.FA;
F_L_N = config.matrices.F_L_d;
FA_Rdd_AF_tx1 = config.matrices.FA_Rdd_AF_tx1;
XpF = config.preamble.XpF_N;
R_dd = config.matrices.R_dd;

off_bins = config.off_bins;
on_bins1 = config.on_bins1;
on_bins_single = config.on_bins_single;

Xp_iT = config.preamble.Xp_uw;

A = config.matrices.A_uw;

varN = 1/(10^(snr/10));

% Channel Equalization
X_hat = zeros(N,nT);
equalizer = zeros(nT*N,nR*N);
equalizer_cov = zeros(nT*N,nR*N);
for n = 1:N
    H = squeeze(H_hat(n,:,:)).';
    equalizer(n:N:end,n:N:end) = (H'/(H*H' + varN*eye(nR) ));    
    X_hat(n,:) = equalizer(n:N:end,n:N:end)*Y_all(n,:).';
    
    equalizer_cov(n:N:end,n:N:end) = equalizer(n:N:end,n:N:end)*H; 
end
x_hat = ifft_u(X_hat);    
x_hat = x_hat(1:N_data,:);


A_equalizer = A'/(A*A' + varN*eye(N_data));

d_hat_nT = zeros(N_data,nT);
for iT = 1:nT    
    d_hat_nT(:,iT) = A_equalizer*x_hat(:,iT);
end

d_hat = d_hat_nT(:);

F_NNd = dftmtx(N)/sqrt(N);
F_NNd = F_NNd(:,1:N_data);



R_dd_hat = kron(eye(nT),A_equalizer*F_NNd')*sparse(equalizer_cov)*kron(eye(nT),F_NNd*A);

d_hat = d_hat./diag(R_dd_hat);

R_dd_tilde = eye(nT*N_data) - R_dd_hat;




    
    

end

