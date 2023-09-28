function [ d_hat, h_hat_out, Rx_matrices ] = do_uw_joint_eq_genie_ce(config, snr, y, h_hat_primary, R_hh_err_primary, noblk)
% Joint Channel Estimation and Equalization for Unique Words
%
% Shahab Ehsanfar, TU Dresden


nT = config.nT;
nR = config.nR;
K = config.payload.p.K;
M = config.payload.p.M;
p = config.payload.p;
N_data = config.payload.N_data; 
Np = config.Np;
numTaps = config.numTaps; % number of channel taps
fd = config.fd; % maximum Doppler shift
Fs = config.Fs; % Sampling frequency
n_H = config.n_H; % sample index where CE is at its best
N = config.N; % N = N_data + Np = N_data + N_preamble/2
B = 8;

n_Hd = Np + N_data/2; 

h_hat_out = cell(noblk,1);
% h_hat_out{1} = h_hat_primary;
h_hat_all = h_hat_primary;

h_hat_in = h_hat_primary;
R_hh_err_in = R_hh_err_primary;


d_hat = cell(noblk,1);
for blk = 1:noblk  
    
    

    % Genie-aided CE
        h_hat_in = config.h_genie{blk};
%         h_hat_all = cat(4, h_hat_all, h_hat_secondary);
        
        R_hh_err_in = zeros(nT*numTaps,nT*numTaps,nR);
        
        h_hat_out{blk} = h_hat_in;
 
    
    
    % CP Emulation at the receiver
    yd = do_emulate_cp(config,h_hat_in,y,blk,'for_data');
    
    % FFT Operation
    Y_all = fft_u(yd);
    Y_all_v = Y_all(:);
    H_hat_in = zeros(N_data,nT,nR);
    for iR = 1:nR
        H_hat_in(:,:,iR) = fft(h_hat_in(:,:,iR),N_data);
    end
    
    % UW MMSE Equalization
    [d_hat{blk}, ~, ~, R_dd_hat, ~] = do_mimo_uw_mmse_equalization(config, snr, H_hat_in, R_hh_err_in, Y_all_v);

        
%     % Interference removal for secondary channel estimation
%     d_hat_mtx = reshape(d_hat{blk},[length(d_hat{blk})/nT nT]);
%     Dd_hat = zeros(K,M,nT);
%     xd_hat = zeros(N_data,nT);
%     for iT = 1:nT
%         Dd_hat (:,:,iT) = reshape(d_hat_mtx(:,iT),[K M]);
%         xd_hat(:,iT) = do_modulate(p,Dd_hat(:,:,iT));
%     end    
%     %yp = do_emulate_cp(config,h_hat_in,y,blk,'for_uw',xd_hat);
        
    
    
   % Rx_matrices.R_dd_tilde{blk} = R_dd_tilde;
    
    Rx_matrices.R_dd_hat{blk} = diag(R_dd_hat);
    

end


end
