function [ d_hat, h_hat_out, Rx_matrices ] = do_uw_joint_ce_eq_circUW(config, snr, y, h_hat_primary, R_hh_err_primary, noblk, wiener_flag)
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

N = config.N; % N = N_data + Np = N_data + N_preamble/2
n_H = config.n_H+(N+Np)/2; % sample index where CE is at its best
B = 8;

n_Hd = Np + N_data/2; 

h_hat_out = cell(noblk,1);
% h_hat_out{1} = h_hat_primary;
h_hat_all = [];%h_hat_primary;

h_hat_in = h_hat_primary;
R_hh_err_in = R_hh_err_primary;


d_hat = cell(noblk,1);
for blk = 1:noblk  
    
    % Take UW from UW-Payload-UW
    yp = y([numTaps:Np end-Np+numTaps:end],:,blk);
    config.ce_blk = blk;
    
    % Secondary Channel Estimation
    h_hat_secondary = zeros(numTaps,nT,nR);
    R_hh_err_2nd = zeros(nT*numTaps,nT*numTaps,nR);
%     if (wiener_flag ~= 2 || blk <= B)
        for iR = 1:nR
            %[h_hat_secondary(:,:,iR), R_hh_err_2nd(:,:,iR)] = do_mimo_uw_channel_estimation( config, snr, yp(:,iR), 'secondary_Slmmse', h_hat_in(:,:,iR),R_hh_err_in(:,:,iR),diag(diag(R_dd_tilde)));
            %[h_hat_secondary(:,:,iR), R_hh_err_2nd(:,:,iR)] = do_mimo_uw_channel_estimation( config, snr, yp(:,iR), 'secondary_lmmse', h_hat_in(:,:,iR),R_hh_err_in(:,:,iR),diag(diag(R_dd_tilde)));
            [h_hat_secondary(:,:,iR), R_hh_err_2nd(:,:,iR)] = do_mimo_uw_channel_estimation( config, snr, yp(:,iR), 'secondary_lmmse_skipL_circUW');
            %[h_hat_secondary(:,:,iR), R_hh_err_2nd(:,:,iR)] = do_mimo_uw_channel_estimation( config, snr, yp(:,iR), 'secondary_slmmse_skipL',h_hat_in(:,:,iR),R_hh_err_in(:,:,iR));
        end
%     elseif blk > B
%         for iR = 1:nR
%             [h_hat_secondary(:,:,iR), R_hh_err_2nd(:,:,iR)] = do_mimo_uw_channel_estimation( config, snr, yp(:,iR), 'secondary_slmmse_skipL',h_hat_in(:,:,iR),R_hh_err_in(:,:,iR));
%         end
%     end
    

    % WIENER FILTERING
    if (wiener_flag == 1 && blk > 1)
        h_all_blks = zeros(numTaps,nT,nR,blk);
        Covariance = zeros(numTaps,numTaps,nT,nR,blk);
        for iR = 1:nR
            for iT = 1:nT
                h_hat_wiener_in = cat(2, squeeze(h_hat_all(:,iT,iR,:)), h_hat_secondary(:,iT,iR));
                R_hh_err_2nd_iT = R_hh_err_2nd((iT-1)*numTaps+1:iT*numTaps,(iT-1)*numTaps+1:iT*numTaps,iR);
                [h_all_blks(:,iT,iR,:), Covariance(:,:,iT,iR,:)] = do_wiener_smoodiction_CircUW(h_hat_wiener_in, R_hh_err_2nd_iT, size(h_hat_wiener_in,2), N, n_H, n_Hd, Fs, fd, numTaps);
            end
        end
        h_hat_in = h_all_blks(:,:,:,end);
        
        h_hat_all = cat(4, h_hat_all, h_hat_secondary);
        %     h_hat_in = h_hat_secondary;
        for iR = 1:nR
            if nT == 4
            R_hh_err_in(:,:,iR) = blkdiag((Covariance(:,:,1,iR,end)),(Covariance(:,:,2,iR,end)), ...
                (Covariance(:,:,3,iR,end)),(Covariance(:,:,4,iR,end)));
            elseif nT == 1
                R_hh_err_in(:,:,iR) = Covariance(:,:,1,iR,end);
            end
        end
        
        h_hat_out{blk} = h_hat_in;
%     elseif (wiener_flag == 2 && blk >= B)
%         h_all_blks = zeros(numPaths,nT,nR,blk+3);
%         Covariance = zeros(numPaths,numPaths,nT,nR,blk+3);
%         for iR = 1:nR
%             for iT = 1:nT
%                 h_hat_wiener_in = cat(2, squeeze(h_hat_all(:,iT,iR,:)), h_hat_secondary(:,iT,iR));
%                 R_hh_err_2nd_iT = R_hh_err_2nd((iT-1)*numPaths+1:iT*numPaths,(iT-1)*numPaths+1:iT*numPaths,iR);
%                 [h_all_blks(:,iT,iR,:), Covariance(:,:,iT,iR,:)] = do_wiener_smoodiction_uw(h_hat_wiener_in, R_hh_err_2nd_iT, size(h_hat_wiener_in,2)+1, N, n_H, Fs, fd, numPaths);
%             end
%         end
%         h_hat_in = h_all_blks(:,:,:,end);
%         
%         h_hat_all = cat(4, h_hat_all, h_all_blks(:,:,:,end-1));
%         %     h_hat_in = h_hat_secondary;
%         for iR = 1:nR
%             R_hh_err_in(:,:,iR) = blkdiag((Covariance(:,:,1,iR,end)),(Covariance(:,:,2,iR,end)), ...
%                 (Covariance(:,:,3,iR,end)),(Covariance(:,:,4,iR,end)));
%         end
%         
%         h_hat_out{blk+1} = h_all_blks(:,:,:,end-1);
    else
        h_hat_in =  h_hat_secondary;
        h_hat_all = cat(4, h_hat_all, h_hat_secondary);
        
        R_hh_err_in = R_hh_err_2nd;
        
        h_hat_out{blk} = h_hat_in;
    end
    
    
    % CP Emulation at the receiver
    yd = do_emulate_cp_circUW(config,h_hat_in,y,blk,'for_data');
    
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
