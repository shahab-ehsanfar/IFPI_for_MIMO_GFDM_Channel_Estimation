function [ H_hat_out, R_hh_err ] = do_mimo_uw_channel_estimation( config, snr, rx_signal, mode, h_hat_in, R_hh_err_in, R_dd_tilde)
% Channel Estimation for Unique-Word transmission scheme
% 


if strcmp(mode, 'primary')
    XpF = config.preamble.XpF;
    FXp_XpF = config.preamble.FXp_XpF;
    R_hh = config.preamble.R_hh;
    nT = config.nT;
    
    sigma_2 = 10^(-snr/10);
    
    % MMSE Channel Estimation
    MMSE_estimator = (sigma_2*inv(R_hh) + FXp_XpF)\XpF';
    
    Yp = fft_u(rx_signal);
    
    h_hat = MMSE_estimator*Yp;
        
    R_hh_err = R_hh - MMSE_estimator*XpF*R_hh;    
       
    % h vector to matrix
    h_hat_mtx = reshape(h_hat,[length(h_hat)/nT nT]);
    
    H_hat_out = h_hat_mtx;
%     % Interpolation
%     H_hat_out = fft(h_hat_mtx,frame_length2);
elseif strcmp(mode, 'primary_toeplitz')
    xp_NpreambleminusL_L = config.preamble.xp_NpreambleminusL_L;
    R_hh = config.preamble.R_hh;
    nT = config.nT;
    numPaths = config.numPaths;
    
    sigma_2 = 10^(-snr/10);
    
    MMSE_estimator = (sigma_2*inv(R_hh) + xp_NpreambleminusL_L'*xp_NpreambleminusL_L)\xp_NpreambleminusL_L';
    
    yp = rx_signal;%rx_signal(numPaths:end);
    
    h_hat = MMSE_estimator*yp;
    
    R_hh_err = R_hh - MMSE_estimator*xp_NpreambleminusL_L*R_hh;    
       
    % h vector to matrix
    h_hat_mtx = reshape(h_hat,[length(h_hat)/nT nT]);
    
    H_hat_out = h_hat_mtx;
    
%     % Interpolation
%     H_hat_out = fft(h_hat_mtx,frame_length2);
elseif strcmp(mode, 'primary_toeplitz_end')
    xp_NpreambleminusL_endL = config.preamble.xp_NpreambleminusL_endL;
    R_hh_endL = config.preamble.R_hh_endL;
    nT = config.nT;
    
    
    sigma_2 = 10^(-snr/10);
    
    MMSE_estimator = (sigma_2*inv(R_hh_endL) + xp_NpreambleminusL_endL'*xp_NpreambleminusL_endL)\xp_NpreambleminusL_endL';
    
    yp = rx_signal;%rx_signal(numPaths:end);
    
    h_hat = MMSE_estimator*yp;
    
    R_hh_err = R_hh_endL - MMSE_estimator*xp_NpreambleminusL_endL*R_hh_endL;
    
    % h vector to matrix
    h_hat_mtx = reshape(h_hat,[length(h_hat)/nT nT]);
    
    H_hat_out = h_hat_mtx;
    
%     % Interpolation
%     H_hat_out = fft(h_hat_mtx,frame_length2);
elseif strcmp(mode, 'secondary_slmmse_skipL')
    xp_NpminusL_L = config.preamble.xp_NpminusL_L;
%     R_hh = config.preamble.R_hh;
    nT = config.nT;
    numPaths = config.numPaths;
    
    sigma_2 = 10^(-snr/10);
    
    yp = rx_signal(numPaths:end) - xp_NpminusL_L*h_hat_in(:);
    
    
%     R_yy = xp_NpminusL_L*R_hh_err_in*xp_NpminusL_L' + sigma_2*eye(size(xp_NpminusL_L,1));
%     R_hy = R_hh_err_in*xp_NpminusL_L';
%     MMSE_estimator = R_hy/R_yy;
%     R_hh_err = R_hh_err_in - MMSE_estimator*R_hy'; 
    
    
    MMSE_estimator = (sigma_2*inv(R_hh_err_in) + xp_NpminusL_L'*xp_NpminusL_L)\xp_NpminusL_L';
    R_hh_err = R_hh_err_in - MMSE_estimator*xp_NpminusL_L*R_hh_err_in;    
    
   
    
    h_hat = h_hat_in(:) + MMSE_estimator*yp;
    
       
       
    % h vector to matrix
    h_hat_mtx = reshape(h_hat,[length(h_hat)/nT nT]);
    
    H_hat_out = h_hat_mtx;

elseif strcmp(mode, 'secondary_lmmse_skipL')
    xp_NpminusL_L = config.preamble.xp_NpminusL_L;
    R_hh = config.preamble.R_hh;
    nT = config.nT;
    numPaths = config.numPaths;
    
    sigma_2 = 10^(-snr/10);
    
    MMSE_estimator = (sigma_2*inv(R_hh)*0+sigma_2*eye(size(R_hh,1)) + xp_NpminusL_L'*xp_NpminusL_L)\xp_NpminusL_L';
    
    yp = rx_signal(numPaths:end);
    
    h_hat = MMSE_estimator*yp;
    
    R_hh_err = R_hh - MMSE_estimator*xp_NpminusL_L*R_hh;    
       
    % h vector to matrix
    h_hat_mtx = reshape(h_hat,[length(h_hat)/nT nT]);
    
    H_hat_out = h_hat_mtx;
elseif strcmp(mode, 'secondary_lmmse_skipL_circUW')
    blk = config.ce_blk;    
    xp_NpminusL_L = config.preamble.xp_NpminusL_L;
    xp_NpminusL_L = circshift(xp_NpminusL_L,[(size(xp_NpminusL_L,1)/2)*(1-mod(blk,2)) 0]);
    R_hh = config.preamble.R_hh;
    nT = config.nT;
%     numPaths = config.numPaths;
    
    sigma_2 = 10^(-snr/10);
    
    MMSE_estimator = (sigma_2*inv(R_hh) + xp_NpminusL_L'*xp_NpminusL_L)\xp_NpminusL_L';
    
    yp = rx_signal;
    
    h_hat = MMSE_estimator*yp;
    
    R_hh_err = R_hh - MMSE_estimator*xp_NpminusL_L*R_hh;    
       
    % h vector to matrix
    h_hat_mtx = reshape(h_hat,[length(h_hat)/nT nT]);
    
    H_hat_out = h_hat_mtx;    
elseif strcmp(mode, 'secondary_lmmse')
    xp_NpminusL_L = config.preamble.xp_NpminusL_L;
    R_hh = config.preamble.R_hh;
    nT = config.nT;
    numPaths = config.numPaths;
    
    sigma_2 = 10^(-snr/10);
    
    MMSE_estimator = (sigma_2*inv(R_hh) + xp_NpminusL_L'*xp_NpminusL_L)\xp_NpminusL_L';
    
    yp = rx_signal;
    
    h_hat = MMSE_estimator*yp;
    
    R_hh_err = R_hh - MMSE_estimator*xp_NpminusL_L*R_hh;    
       
    % h vector to matrix
    h_hat_mtx = reshape(h_hat,[length(h_hat)/nT nT]);
    
    H_hat_out = h_hat_mtx;
    
elseif strcmp(mode, 'secondary_lmmse_approx')
    XpF = config.preamble.XpF;
    FXp_XpF = config.preamble.FXp_XpF;
    R_hh = config.preamble.R_hh;
    nT = config.nT;    
        
    sigma_2 = 10^(-snr/10);
    
    % MMSE Channel Estimation
    MMSE_estimator = (sigma_2*inv(R_hh) + FXp_XpF)\XpF';
    
    Yp = fft_u(rx_signal);
    
    h_hat = MMSE_estimator*Yp;
    
    
    R_hh_err = R_hh - MMSE_estimator*XpF*R_hh;    
       
    % h vector to matrix
    h_hat_mtx = reshape(h_hat,[length(h_hat)/nT nT]);
    
    H_hat_out = h_hat_mtx;
        
    
    
elseif strcmp(mode, 'secondary_lmmse_old')    
    XpF = config.preamble.XpF;    
    R_hh = config.preamble.R_hh;
    nT = config.nT;      
    F_Np_F_NdNp = config.matrices.F_Np_F_NdNp;
    FA_Rdd_AF_tx1 = config.matrices.FA_Rdd_AF_tx1;
    FA = config.matrices.FA;
    F_L_d = config.matrices.F_L_d;
    N_data = config.payload.N_data;
    numPaths = config.numPaths;
    Np = config.Np;
    
    Yp = fft_u(rx_signal);
    
    sigma_2 = 10^(-snr/10);
    
    H_hat_primary = zeros(N_data,nT);
    R_XX_hat = zeros(N_data,N_data,nT);
    for iT= 1:nT
    H_hat_primary(:,iT) = fft(h_hat_in(:,iT),N_data);
    R_XX_hat(:,:,iT) = FA*R_dd_tilde((iT-1)*N_data+1:iT*N_data,(iT-1)*N_data+1:iT*N_data)*FA';    
    end
    
    mat2_inner = 0;
    mat3_inner = 0;
    for iT = 1:nT              
         R_HH_hat = F_L_d*R_hh_err_in((iT-1)*numPaths+1:iT*numPaths,(iT-1)*numPaths+1:iT*numPaths)*F_L_d';
         mat2_inner = mat2_inner + (R_HH_hat.*FA_Rdd_AF_tx1);         
         mat3_inner = mat3_inner + repmat(H_hat_primary(:,iT), [1 N_data]).*R_XX_hat(:,:,iT).*repmat(H_hat_primary(:,iT)', [N_data 1]);
    end
    
    mat1 = XpF*R_hh*XpF';
    mat2 = F_Np_F_NdNp*mat2_inner*F_Np_F_NdNp';    
    mat3 = F_Np_F_NdNp*mat3_inner*F_Np_F_NdNp';
    mat4 = 2*sigma_2*eye(Np);
    
    
    R_yy = mat1+mat2+mat3+mat4;
    R_hy = R_hh*XpF';
    
    MMSE_estimator = R_hy/R_yy;
    h_hat = MMSE_estimator*Yp;
        
    
    
    R_hh_err = R_hh - MMSE_estimator*XpF*R_hh;    
       
    % h vector to matrix
    h_hat_mtx = reshape(h_hat,[length(h_hat)/nT nT]);
    
    H_hat_out = h_hat_mtx;
    
elseif strcmp(mode, 'secondary_Slmmse')    
    XpF = config.preamble.XpF;        
    nT = config.nT;      
    F_Np_F_NdNp = config.matrices.F_Np_F_NdNp;
    FA_Rdd_AF_tx1 = config.matrices.FA_Rdd_AF_tx1;
    FA = config.matrices.FA;
    F_L_d = config.matrices.F_L_d;
    N_data = config.payload.N_data;
    numPaths = config.numPaths;
    Np = config.Np;
    R_hh_err0 = R_hh_err_in;
    
    Yp = fft_u(rx_signal);
    
    Yp_hat = XpF*h_hat_in(:);
    
    Yp_tilde = Yp - Yp_hat;
    
    
    sigma_2 = 10^(-snr/10);
    
    H_hat_primary = zeros(N_data,nT);
    R_XX_hat = zeros(N_data,N_data,nT);
    for iT= 1:nT
    H_hat_primary(:,iT) = fft(h_hat_in(:,iT),N_data);
    R_XX_hat(:,:,iT) = FA*R_dd_tilde((iT-1)*N_data+1:iT*N_data,(iT-1)*N_data+1:iT*N_data)*FA';    
    end
    
    mat2_inner = 0;
    mat3_inner = 0;
    for iT = 1:nT              
         R_HH_hat = F_L_d*R_hh_err0((iT-1)*numPaths+1:iT*numPaths,(iT-1)*numPaths+1:iT*numPaths)*F_L_d';
         mat2_inner = mat2_inner + (R_HH_hat.*FA_Rdd_AF_tx1);         
         mat3_inner = mat3_inner + repmat(H_hat_primary(:,iT), [1 N_data]).*R_XX_hat(:,:,iT).*repmat(H_hat_primary(:,iT)', [N_data 1]);
    end
    
    mat1 = XpF*R_hh_err0*XpF';
    mat2 = F_Np_F_NdNp*mat2_inner*F_Np_F_NdNp';    
    mat3 = F_Np_F_NdNp*mat3_inner*F_Np_F_NdNp';
    mat4 = 2*sigma_2*eye(Np);
    
    
    R_yy = mat1+mat2+mat3+mat4;
    R_hy = R_hh_err0*XpF';
    
    MMSE_estimator = R_hy/R_yy;
    h_hat_tilde = MMSE_estimator*Yp_tilde;
        
    h_hat = h_hat_in(:) + h_hat_tilde;
    
    R_hh_err = R_hh_err0 - MMSE_estimator*XpF*R_hh_err0;    
       
    % h vector to matrix
    h_hat_mtx = reshape(h_hat,[length(h_hat)/nT nT]);
    
    H_hat_out = h_hat_mtx;

elseif strcmp(mode, 'other')
    
    nT = config.nT;
    nR = config.nR;
    N_blk = size(H_hat_in,1);
    Xp_uw = config.preamble.Xp_uw;   
    
    
    Yp_iR = zeros(N_blk,nR);
    for iR = 1:nR
       Yp_iR(:,iR) = fft_u(rx_signal(:,iR)); 
    end
 
%     % TO BE DELETED
%     % Reference for equalization
%     H_all = [];    
%     for iR = 1:nR       
%         Yp_iR(:,iR) = fft_u(rx_signal(:,iR));
%         
%         H_row = [];
%         for iT = 1:nT
%            H_row = [H_row diag(H_hat_in(:,iT,iR))];           
%         end
%         H_all = [H_all; H_row];   
%     end    
%     Xd_hat0 = (H_all'*H_all)\H_all'*Yp_iR(:);
%     Xd_hat_mtx = reshape(Xd_hat0,[length(Xd_hat0)/nT nT]);
%     xd_hat_ifft = ifft(Xd_hat_mtx); figure; plot(abs(xd_hat_ifft))
    
    

    % Channel Equalization for interference removal
    X_hat = zeros(N_blk,nT);
    for n = 1:N_blk
       H_n = squeeze(H_hat_in(n,:,:)).';       
%        cond_H(n) = cond(H_n'*H_n);
%        if cond_H(n) > 30
%            H_n_H = (H_n'*H_n + eye(nR) );
%        else
%            H_n_H = (H_n'*H_n);
%        end
       H_n_H = (H_n'*H_n + 10^(-snr/10));
       equalizer = H_n_H\H_n';
       X_hat_biased = equalizer*(Yp_iR(n,:).');
       X_hat(n,:) = X_hat_biased./diag(H_n*equalizer);
       
       
%        H_n_H = (H_n*H_n' + 10^(-snr/10));
%        X_hat(n,:) = (H_n'/H_n_H)*(Yp_iR(n,:).');
       
       cond_new(n) = cond(H_n_H);
%        X_hat(n,:) = H_n\(Yp_iR(n,:).');
    end
    figure; plot(abs(cond_new))
    figure; plot(abs(X_hat(:)))
    figure; plot(abs(ifft(X_hat)))
    
    interference = X_hat - Xp_uw;
    figure; plot(abs(ifft_u(interference)));
    []
    
    
    
    
    sigma_2 = 10^(-snr/10);
    N_data = config.payload.N_data;
    N_preamble = config.preamble.N_preamble;
    R_d = blkdiag(eye(N_data),zeros(N_preamble/2));
    F_big = dftmtx(N_data+N_preamble/2)/sqrt(N_data+N_preamble/2);    
    FRdF = F_big*R_d*F_big';
    
end





end

