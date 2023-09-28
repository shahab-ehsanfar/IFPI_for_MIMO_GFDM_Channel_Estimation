function [ vH_lmmse_out, d_hat ] = do_CE_MMSE_PIC( p, conf, pilots_conf, no_PIC_iterations, snr, vH_lmmse_in, R_hh_err, Dd1_indices, ycp, Y_all_v, dhat_genie)
% Parallel Interference Cancellation for MIMO GFDM Channel Estimation
%
% TU-Dresden
% Shahab Ehsanfar
%

nT = conf.nT;
nR = conf.nR;
N = conf.N;
N_tot = conf.N_tot;
N_2blk = conf.N_2blk;
Ncp = conf.Ncp;
Fbig_cp = conf.Fbig_cp;
F_Lcp = conf.F_Lcp;
P = conf.P;
numPaths = conf.numPaths;
Delta_k = conf.Delta_k;
FA = conf.FA;
F_L_N = conf.F_L_N;
JF = conf.JF;
noblk = conf.noblk;
FIA = conf.FIA;

off_bins = pilots_conf.off_bins;
on_bins1 = pilots_conf.on_bins1;
on_bins_single = pilots_conf.on_bins_single;
FA_Rdd_AF_tx1 = pilots_conf.FA_Rdd_AF_tx1;
XpF = pilots_conf.XpF;
R_dd = pilots_conf.R_dd;
XpF_cp = pilots_conf.XpF_cp;
Xp_iT = pilots_conf.Xp_iT;


for its=1:no_PIC_iterations
    
    
    wH_lmmse_in = reshape(vH_lmmse_in(:,its), [N nT*nR]);    
   % wH_lmmse_in = reshape(conf.vH_big, [N nT*nR]);    
    H_long = fft(ifft(wH_lmmse_in),(N_2blk));   
       
    
    % HA = H_hat*kron(eye(nT),F*A);
    HA = [];
    for iR = 1:nR
        HA_iR = [];
        for iT = 1:nT
            HA_iR = [HA_iR (repmat(wH_lmmse_in(:,iT+(iR-1)*nR),[1 N]).*FA)];
        end
        HA = [HA; HA_iR];
    end
    HA(off_bins,:) = [];
    
    
    R_hhe_i = zeros(numPaths*nT,numPaths*nT,nR);
    R_HHe_i = zeros(N,N,nT);
    %    R_HHe_i_cpbig = zeros(N_2blk,N_2blk,nT,nR); R_HHe_i_cpbig_cor12 = zeros(N_2blk,N_2blk,nT,nR);R_HHe_i_cpbig_cor21 = zeros(N_2blk,N_2blk,nT,nR);
    R_psipsi_i = zeros(N,N,nR);
    R_psidp_i = zeros(N,N,nR);
    for iR = 0:nR-1
        indx = (numPaths*nT*iR);
        R_hhe_i(:,:,iR+1) = R_hh_err(indx+1:indx+numPaths*nT,indx+1:indx+numPaths*nT);
        % R_HHe_i(:,:,iR+1) =  (N)*F_L_w*R_hhe_i(:,:,iR+1)*F_L_w';
        for iT = 0:nT-1
            R_HHe_i(:,:,iT+1) = N*F_L_N*R_hhe_i(iT*numPaths+1:numPaths+iT*numPaths,iT*numPaths+1:numPaths+iT*numPaths,iR+1)*F_L_N';
            
            R_psipsi_i(:,:,iT+1) = FA_Rdd_AF_tx1 .* R_HHe_i(:,:,iT+1);
            
            %                     R_HHe_i_cpbig(:,:,iT+1,iR+1) = N_2blk*F_Lcpblk*R_hhe_i(iT*numPaths+1:numPaths+iT*numPaths,iT*numPaths+1:numPaths+iT*numPaths,iR+1)*F_Lcpblk';
            %                     if iT == 0
            %                         R_HHe_i_cpbig_cor12(:,:,iT+1,iR+1) = N_2blk*F_Lcpblk*R_hhe_i(iT*numPaths+1:numPaths+iT*numPaths,(iT+1)*numPaths+1:numPaths+(iT+1)*numPaths,iR+1)*F_Lcpblk';
            %                         R_HHe_i_cpbig_cor21(:,:,iT+1,iR+1) = N_2blk*F_Lcpblk*R_hhe_i((iT+1)*numPaths+1:numPaths+(iT+1)*numPaths,iT*numPaths+1:numPaths+iT*numPaths,iR+1)*F_Lcpblk';
            %                     end
        end
        R_psidp_i(:,:,iR+1) = N*XpF*R_hhe_i(:,:,iR+1)*XpF';
        
        % R_psidp_i(:,:,iR+1) = (FIA*(dpdp)*FIA') .* R_HHe_i(:,:,iR+1);
        % R_psipsi_i(:,:,iR+1) = FA_Rdd_AF_tx1 .* R_HHe_i(:,:,iR+1) % FIARddAIF .* R_HHe_i(:,:,iR+1);
    end
    if nR == 2
        R_psidp_tot = blkdiag(R_psidp_i(on_bins_single,on_bins_single,1),R_psidp_i(on_bins_single,on_bins_single,2));
    elseif nR == 4
        R_psidp_tot = blkdiag(R_psidp_i(on_bins_single,on_bins_single,1),R_psidp_i(on_bins_single,on_bins_single,2),...
            R_psidp_i(on_bins_single,on_bins_single,3),R_psidp_i(on_bins_single,on_bins_single,4));
    end
    R_psipsi_tot = kron(eye(nR),sum(R_psipsi_i(on_bins_single,on_bins_single,:),3));
    
    % R_HH_err =  diag([zeros(66,1); ones(300,1); zeros(66,1)])*((N+2*Ncp))*Fbig_Lcpblk*R_hh_err*Fbig_Lcpblk'*diag([zeros(66,1); ones(300,1); zeros(66,1)]);
    % R_HH_err =  ((N_2blk))*Fbig_Lcpblk_nTnR*R_hh_err*Fbig_Lcpblk_nTnR';
    
    
    
    %          mat1 = 0 ; %dJF*( R_psidp_tot)*dJF';
    %          mat2 = JH_hatIA*R_dd*JH_hatIA';
    %          %mat3 = dJF*( (FIA*R_dd*FIA') .* R_HH_err  )*dJF';
    %          mat3 = 0; %dJF*( R_psipsi_tot)*dJF';
    %
    varN = 1/(10^(snr/10));
    %
    %          R_yy = mat1 + mat2 + mat3 + varN*eye(N*nR);
    %          R_dy = R_dd*(JH_hatIA)';
    %         % R_dy = (R_dd)*(JH_hatIA)';
    %
    %          equalizer = R_dy/R_yy;
    %       %   d_hat_better = equalizer*y(:);
    
    % R_dy = R_dd*HA';
    R_dy = repmat(diag(R_dd),[1 length(on_bins1)]) .* HA';
    
    % mat1 = HA*R_dd*HA';
    mat1 = HA*R_dy;
    mat2 = R_psipsi_tot;
    mat3 = R_psidp_tot;
    
    R_yy = mat1 + mat2 + mat3 + varN*eye(length(on_bins1));
    
    % R_yy = mat1 + varN*eye(length(on_bins1));
    
    % HXp = H_hat*Xp_iT(:);      %*kron(eye(nT),F*A)*dp;
    HXp_mat = zeros(N,nR);
    for iR = 1:nR
        HXp_iR = 0;
        for iT = 1:nT
            HXp_iR = HXp_iR + wH_lmmse_in(:,iT+(iR-1)*nR) .* Xp_iT(:,iT);
        end
        HXp_mat(:,iR) = HXp_iR;
    end
    HXp = HXp_mat(:);
    
    % equalizer = (HA'*HA+ varN*eye((N- length(pilot_indx) )*nR) )\(HA'); %*(Y_all(:) - H_hat*kron(eye(nT),F*A)*dp )
    
    %dp = Dp0(:);
    equalizer = R_dy/R_yy;
    d_hat = equalizer*(Y_all_v(on_bins1) - HXp(on_bins1) );
    
    d_hat_test = do_CEPIC_equalization( conf, pilots_conf, its, snr, vH_lmmse_in, R_hh_err, Y_all_v);
    
 %   d_hat_save(:,its)= d_hat;
    %         mean(abs( Dd(:) - d_hat ).^2 )
    
    %          equalizer2 = R_dy/(mat1 +varN*eye(N*nR));
    %          d_hat2 = equalizer2*(Y_all(:) - HXp );
    %
    %
    %          mean(abs( Dd(:) - d_hat2 ).^2 )
    
    
    %      d_hat_better(1:N) = 1;
%     if nT == 2
%         d_hat_final = [zeros(N,1); d_hat(1:N); zeros(N,1); d_hat(N+1:2*N)];
%         
%     elseif nT == 4
%         d_hat_final = [zeros(N,1); d_hat(1:N); zeros(N,1); d_hat(N+1:2*N); ...
%             zeros(N,1); d_hat(2*N+1:3*N); zeros(N,1); d_hat(3*N+1:4*N)];
%     end    
    
    
    JFHhat = [];
    for iR = 1:nR
        JFHhat_iR = [];
        for iT = 1:nT
            JFHhat_iR = [JFHhat_iR (JF .* repmat(H_long(:,iT+(iR-1)*nR).',[N+Ncp 1]))];
        end
        JFHhat = [JFHhat; JFHhat_iR];
    end    
    
%     JH_hatIA = [];
%     for iT = 0:nT-1
%         JH_hatIA = [JH_hatIA JFHhat(:,1+(iT*N_2blk):N_2blk+(iT*N_2blk))*FIA];
%     end
    
    JH_hatIA_3d = zeros((N+Ncp)*nR,N,nT);
    for iT = 0:nT-1 
       
        JH_hatIA_3d(:,:,iT+1) = JFHhat(:,1+(iT*N_2blk):N_2blk+(iT*N_2blk) )*FIA;
        
    end
    JH_hatIA = reshape(JH_hatIA_3d,[(N+Ncp)*nR N*nT]);
    
    if its > 1
        y2 = ycp(:) - JH_hatIA*dhat_genie;
        R_dd_tilde = zeros(size(R_dd));
    else
        y2 = ycp(:) - JH_hatIA*d_hat; %_final;
        R_dd_tilde = R_dd - equalizer*R_dy';
    end
    
    
    
    %             R_dd_tilde = R_dd + dpdp - equalizer*R_dy';
    
    
    
    y2 = reshape(y2, [N_tot nR]);
    
    
    dR_dd_tilde = diag(R_dd_tilde);
    if nT*nR == 4
        sigma_dhat1 = dR_dd_tilde(1:N);
        sigma_dhat2 = dR_dd_tilde(N+1:2*N);
        Dd_tilde(:,:,1) = reshape(sigma_dhat1,[p.K p.M]);
        Dd_tilde(:,:,2) = reshape(sigma_dhat2,[p.K p.M]);
        
        %             R_dd_tilde_big11 = blkdiag(eye(N),R_dd_tilde(1:N,1:N));
        %             R_dd_tilde_big22 = blkdiag(eye(N),R_dd_tilde(N+1:2*N,N+1:2*N));
        %             R_dd_tilde_big12 = blkdiag(eye(N),R_dd_tilde(1:N,N+1:2*N));
        %             R_dd_tilde_big21 = blkdiag(eye(N),R_dd_tilde(N+1:2*N,1:N));
        %
        %             % test = [R_dd_tilde(1:N,1:N) R_dd_tilde(1:N,N+1:2*N); R_dd_tilde(N+1:2*N,1:N) R_dd_tilde(N+1:2*N,N+1:2*N)];
        %             R_dd_tilde_big = [R_dd_tilde_big11 R_dd_tilde_big12; R_dd_tilde_big21 R_dd_tilde_big22];
        %
        %
        %             FIA_half = F_cpblk*blkdiag(A(end-Ncp+1:end,:),I_cp*A);
        %             R_X1X1 = FIA_half*R_dd_tilde_big11*FIA_half';
        %             R_X1X2 = FIA_half*R_dd_tilde_big12*FIA_half';
        %             R_X2X1 = FIA_half*R_dd_tilde_big21*FIA_half';
        %             R_X2X2 = FIA_half*R_dd_tilde_big22*FIA_half';
        %
        
        
    elseif nT*nR == 16
        sigma_dhat1 = dR_dd_tilde(1:N);
        sigma_dhat2 = dR_dd_tilde(N+1:2*N);
        sigma_dhat3 = dR_dd_tilde(2*N+1:3*N);
        sigma_dhat4 = dR_dd_tilde(3*N+1:4*N);
        Dd_tilde(:,:,1) = reshape(sigma_dhat1,[p.K p.M]);
        Dd_tilde(:,:,2) = reshape(sigma_dhat2,[p.K p.M]);
        Dd_tilde(:,:,3) = reshape(sigma_dhat3,[p.K p.M]);
        Dd_tilde(:,:,4) = reshape(sigma_dhat4,[p.K p.M]);
    end
    
    %     R_PsiPsi_new1 = JH_hatIA*R_dd_tilde_big*JH_hatIA';
    
    %     R_PsiPsi_new2 = dJF*kron(eye(nR),(R_X1X1 .* R_HHe_i_cpbig(:,:,1,1)) + (R_X2X2 .* R_HHe_i_cpbig(:,:,2,1)) + ( R_X1X2 .* R_HHe_i_cpbig_cor21(:,:,1,1) ) + (R_X2X1 .* R_HHe_i_cpbig_cor12(:,:,1,1) )  )*dJF';
    
    %     R_PsiPsi_new_all = R_PsiPsi_new1 + R_PsiPsi_new2;
    
    R_PsiPsi_new = calc_Rpsipsi('R_PsiPsi_CP_Dhat', nT, nR, conf, Dd_tilde, Dd1_indices);
    
    [vH_lmmse_in(:,its+1), R_hh_err] = do_Chan_Estimation_mimo('LMMSE_CP', y2, snr, N_2blk, nT, nR, Fbig_cp, F_Lcp, 0, XpF_cp, P, numPaths, 0, Delta_k, R_PsiPsi_new, Ncp);
    
    %      test = do_Chan_Estimation_mimo('LMMSE_CP', y2, snr(si), N_2blk, nT, nR, Fbig_cp, F_Lcp, 0, XpF_cp, P, numPaths, Dd, Delta_k, R_PsiPsi_new_all,Ncp);
    
    %      mean(abs(vH_lmmse_in(:,sth+1) - test).^2)
    
    %     mean(abs(vH_big - test).^2)
end

vH_lmmse_out = vH_lmmse_in;


end

