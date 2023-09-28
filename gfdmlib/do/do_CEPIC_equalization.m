function [ d_hat , R_dd_tilde, equalizer, R_dd_hat, HA] = do_CEPIC_equalization( conf, pilots_conf, snr, wH_lmmse_in, R_hh_err, Y_all_v)
% Parallel Interference Cancellation for MIMO GFDM Channel Estimation
%
% TU-Dresden
% Shahab Ehsanfar
%

nT = conf.nT;
nR = conf.nR;
N = conf.N;
numPaths = conf.numPaths;
FA = conf.FA;
F_L_N = conf.F_L_N;


off_bins = pilots_conf.off_bins;
on_bins1 = pilots_conf.on_bins1;
on_bins_single = pilots_conf.on_bins_single;
FA_Rdd_AF_tx1 = pilots_conf.FA_Rdd_AF_tx1;
XpF = pilots_conf.XpF;

if pilots_conf.blk == 1
    R_dd = pilots_conf.R_dd_blk1;
else
    R_dd = pilots_conf.R_dd;
end

Xp_iT = pilots_conf.Xp_iT;

       
    
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
    HA(abs(HA)<1e-7) = 0;
    HA = sparse(HA);
    
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
    if nR == 1
        R_psidp_tot = R_psidp_i(on_bins_single,on_bins_single);
    elseif nR == 2
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
    mat3 = 0;%R_psidp_tot;
    
    R_yy = mat1 + mat2 + mat3 + varN*eye(length(on_bins1));
    R_yy(abs(R_yy)<1e-7) = 0;
    
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
%     tic
%     equalizer2 = R_dy/R_yy;
%     toc
%     tic
%     R_yy(abs(R_yy)<1e-10) = 0;
%     equalizer = R_dy/R_yy;
%     toc

    
    rows = zeros(1,length(on_bins1));
    for iR = 0:nR-1
       rows(iR+1:nR:length(on_bins1)) = (iR*length(on_bins1)/nR+1):(iR+1)*length(on_bins1)/nR;
    end
    % tic
    %spparms('spumoni',2)    
    spparms('bandden',0) % Force Matlab to use banded solver
    equalizer(:,rows) = sparse(R_dy(:,rows))/sparse(R_yy(rows,rows));
    % toc
    
    if pilots_conf.blk == 1
        d_hat = equalizer*(Y_all_v(on_bins1));
    else
        d_hat = sparse(equalizer)*(Y_all_v(on_bins1) - HXp(on_bins1) );
    end
    R_dd_hat = sparse(equalizer(:,rows))*sparse(R_dy(:,rows))';
    R_dd_tilde = sparse(R_dd) - R_dd_hat;
    %R_dd_tilde = R_dd - equalizer*R_dy';
    
    diagR_dd_hat = diag(R_dd_hat);
    diagR_dd_hat(abs(diagR_dd_hat)<1e-4) = 1;
    
    d_hat = d_hat./diagR_dd_hat;
    
    R_dd_tilde = full(R_dd_tilde);
    
    
    
    
%      rows = zeros(1,length(on_bins1));
%      rows(1:2:length(on_bins1)) = 1:length(on_bins1)/2;
%      rows(2:2:length(on_bins1)) = (length(on_bins1)/2+1):length(on_bins1);
%     equalizer_new(:,rows) = R_dy(:,rows)/R_yy(rows,rows);

%     R_yy(abs(R_yy)<1e-10) = 0;
%     figure; spy(abs(R_yy))
%     figure; spy(abs(R_yy(rows,rows)))

end

