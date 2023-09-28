function [ vH_lmmse_out, d_hat, Heq ] = do_CE_MMSE_PIC(est_type, p, conf, pilots_conf, no_PIC_iterations, snr, vH_lmmse_in, R_hh_err, Dd1_indices, ycp, Y_all_v, dhat_genie, vH_genie)
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
if strcmp(est_type,'LMMSE_CP')
    path_delays = conf.path_delays_cpblk;
else
    path_delays = conf.path_delays;
end

P = conf.P;
numPaths = conf.numPaths;
Delta_k = conf.Delta_k;

JF = conf.JF;

FIA = conf.FIA;
FIcpA= conf.FIcpA;

R_dd = pilots_conf.R_dd;


%actual_snr = snr;
%if (strcmp(est_type,'LMMSE_Zm') || strcmp(est_type,'LMMSE_CPZm'))
%    if actual_snr >= 25
%       snr = 15;
%    end
%end
snr_in = snr;
%if length(snr_in) == 2
%   snr = snr(1);
%end

if strcmp(est_type,'LMMSE_CPZm')
	zeta = 1;
else
	zeta = 0;
end

for its=1:no_PIC_iterations+2
    
    if (length(snr_in) == 2 && its <= no_PIC_iterations)
        if (snr_in(2) > snr_in(1))
		snr = snr_in(1) + (its-1)*((snr_in(2)-snr_in(1))/(no_PIC_iterations+zeta)); 
        else
		snr = snr_in(2);
        end
    end

    if its < no_PIC_iterations+1
        wH_lmmse_in = reshape(vH_lmmse_in(:,its), [N nT*nR]);    
   % wH_lmmse_in = reshape(conf.vH_big, [N nT*nR]);    
    elseif its == no_PIC_iterations+1
        wH_lmmse_in = reshape(vH_lmmse_in(:,its-1), [N nT*nR]);    
    else
        wH_lmmse_in = reshape(vH_genie, [N nT*nR]);   
    end
      
   
    if its < no_PIC_iterations+1
        % Second block
        pilots_conf.blk = 2;
        [d_hat, R_dd_tilde, ~, ~, ~] = do_CEPIC_equalization( conf, pilots_conf, snr, wH_lmmse_in, R_hh_err, Y_all_v(:,2));
      %  if strcmp(est_type,'LMMSE_CPZm')
            % First block
            pilots_conf.blk = 1;
            [d_hat1, R_dd_tilde1, ~, ~, ~] = do_CEPIC_equalization( conf, pilots_conf, snr, wH_lmmse_in, R_hh_err, Y_all_v(:,1));
      %  end
    end
    
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
    
    if (strcmp(est_type,'LMMSE_CP') || strcmp(est_type,'LMMSE_CPZm') )
        H_long = fft(ifft(wH_lmmse_in),(N_2blk));
        JFHhat = [];
        for iR = 1:nR
            JFHhat_iR = [];
            for iT = 1:nT
                JFHhat_iR = [JFHhat_iR (JF .* repmat(H_long(:,iT+(iR-1)*nR).',[N+Ncp 1]))];
            end
            JFHhat = [JFHhat; JFHhat_iR];
        end
        
       % if strcmp(est_type,'LMMSE_CPZm')
            JH_hatIA_3d = zeros((N+Ncp)*nR,2*N,nT);
       % else
       %     JH_hatIA_3d = zeros((N+Ncp)*nR,N,nT);
       % end
        for iT = 0:nT-1
           % if strcmp(est_type,'LMMSE_CPZm')
                JH_hatIA_3d(:,:,iT+1) = JFHhat(:,1+(iT*N_2blk):N_2blk+(iT*N_2blk) )*FIcpA;
           % else
           %     JH_hatIA_3d(:,:,iT+1) = JFHhat(:,1+(iT*N_2blk):N_2blk+(iT*N_2blk) )*FIA;
           % end
        end
       % if strcmp(est_type,'LMMSE_CPZm')
            JH_hatIA = reshape(JH_hatIA_3d,[(N+Ncp)*nR 2*N*nT]);
       % else
       %     JH_hatIA = reshape(JH_hatIA_3d,[(N+Ncp)*nR N*nT]);
       % end
    else
        H_long = wH_lmmse_in; 
        JFHhat = [];
        F = dftmtx(N)/sqrt(N);
        for iR = 1:nR
            JFHhat_iR = [];
            for iT = 1:nT
                JFHhat_iR = [JFHhat_iR (F' .* repmat(H_long(:,iT+(iR-1)*nR).',[N 1]))];
            end
            JFHhat = [JFHhat; JFHhat_iR];
        end
        
        JH_hatIA_3d = zeros((N)*nR,N,nT);
        for iT = 0:nT-1
            
            JH_hatIA_3d(:,:,iT+1) = JFHhat(:,1+(iT*N):N+(iT*N) )*conf.FA;
            
        end
        JH_hatIA = reshape(JH_hatIA_3d,[(N)*nR N*nT]);
    end
    
%     JH_hatIA = [];
%     for iT = 0:nT-1
%         JH_hatIA = [JH_hatIA JFHhat(:,1+(iT*N_2blk):N_2blk+(iT*N_2blk))*FIA];
%     end
    
    
    
    
    if its >= no_PIC_iterations+1        
        y2 = ycp(:) - JH_hatIA*dhat_genie;
        R_dd_tilde = zeros(size(R_dd));
    else
        if ( strcmp(est_type,'LMMSE_CPZm') || strcmp(est_type,'LMMSE_CP') )
            D_hat = reshape(d_hat,[N nT]);
            D_hat1 = reshape(d_hat1,[N nT]);
            D_hat_2blk = [D_hat1; D_hat];
            d_hat_cancel = D_hat_2blk(:);
        else
            d_hat_cancel = d_hat;
        end
        y2 = ycp(:) - JH_hatIA*d_hat_cancel; %_final;    JH_hatIA*[zeros(N,1); d_hat(1:N); zeros(N,1); d_hat(N+1:end)]
        %R_dd_tilde = R_dd - equalizer*R_dy';
    end
    
    
    
    %             R_dd_tilde = R_dd + dpdp - equalizer*R_dy';
    
    
    
    
    
    if strcmp(est_type,'LMMSE_CP')
    dR_dd_tilde = diag(R_dd_tilde);
    dR_dd_tilde1 = diag(R_dd_tilde1);
    if nT == 1
    %    sigma_dhat1 = dR_dd_tilde(1:N);
    %    Dd_tilde(:,:,1) = reshape(sigma_dhat1,[p.K p.M]);
        F_L = conf.F_L;
        JFHFIcpA = JFHhat*kron(eye(nR),FIcpA);
                     
     %   R_PsiPsi_new1 = JFHFIcpA*blkdiag(diag(Dd1_indices),diag(diag(R_dd_tilde)))*JFHFIcpA';
        R_PsiPsi_new1 = JFHFIcpA*blkdiag(diag(diag(R_dd_tilde1)),diag(diag(R_dd_tilde)))*JFHFIcpA';
        
        R_HHe_i = zeros(N+2*Ncp,N+2*Ncp,nT);
        R_hhe_i(:,:) = R_hh_err(1:numPaths*nT,1:numPaths*nT);
        iT = 0;
        R_HHe_i(:,:,iT+1) = (N+2*Ncp)*F_L*R_hhe_i(iT*numPaths+1:numPaths+iT*numPaths,iT*numPaths+1:numPaths+iT*numPaths)*F_L';
        
        
        R_X1X1 = FIcpA*blkdiag(diag(Dd1_indices),diag(diag(R_dd(1:N,1:N))))*FIcpA';

        
    elseif nT == 2 
    %    sigma_dhat1 = dR_dd_tilde(1:N);
    %    sigma_dhat2 = dR_dd_tilde(N+1:2*N);
    %    Dd_tilde(:,:,1) = reshape(sigma_dhat1,[p.K p.M]);
    %    Dd_tilde(:,:,2) = reshape(sigma_dhat2,[p.K p.M]);
        
        F_L = conf.F_L;
        
        JFHFIcpA = JFHhat*kron(eye(nR),FIcpA);
        
        %R_PsiPsi_new1 = JFHFIcpA*blkdiag(diag(diag(R_dd_tilde1(1:N,1:N))),diag(diag(R_dd_tilde(1:N,1:N))),diag(diag(R_dd_tilde1(1:N,1:N))),diag(diag(R_dd_tilde(N+1:end,N+1:end))))*JFHFIcpA';
        R_PsiPsi_new1 = (JFHFIcpA.*repmat([diag(R_dd_tilde1(1:N,1:N));diag(R_dd_tilde(1:N,1:N)); ...
                            diag(R_dd_tilde1(N+1:end,N+1:end));diag(R_dd_tilde(N+1:end,N+1:end))].',[nR*N_tot 1]))*JFHFIcpA';
        
        
        R_HHe_i = zeros(N+2*Ncp,N+2*Ncp,nT);
        R_hhe_i(:,:) = R_hh_err(1:numPaths*nT,1:numPaths*nT);
        iT = 0;
        R_HHe_i(:,:,iT+1) = (N+2*Ncp)*F_L*R_hhe_i(iT*numPaths+1:numPaths+iT*numPaths,iT*numPaths+1:numPaths+iT*numPaths)*F_L';
        iT = 1;
        R_HHe_i(:,:,iT+1) = (N+2*Ncp)*F_L*R_hhe_i(iT*numPaths+1:numPaths+iT*numPaths,iT*numPaths+1:numPaths+iT*numPaths)*F_L';
        
        R_X1X1 = FIcpA*blkdiag(diag(Dd1_indices),diag(diag(R_dd(1:N,1:N))))*FIcpA';
                     
                     
        %             R_X1X2 = FIA_half*R_dd_tilde_big12*FIA_half';
        %             R_X2X1 = FIA_half*R_dd_tilde_big21*FIA_half';
        %             R_X2X2 = FIA_half*R_dd_tilde_big22*FIA_half';
        %
        
        
    elseif nT == 4
        F_L = conf.F_L;
        JFHFIcpA = JFHhat*kron(eye(nR),FIcpA);
       
        Dd_tilde_2blk = zeros(2*N,2*N,nT);
        for iT=1:nT
            Dd_tilde_2blk(:,:,iT) = diag([dR_dd_tilde1((iT-1)*N+1: iT*N); dR_dd_tilde((iT-1)*N+1: iT*N) ] );
        end
        Dd_all_nT = get_blkdiag(Dd_tilde_2blk);
        R_PsiPsi_new1 = JFHFIcpA*Dd_all_nT*JFHFIcpA';
        
        R_HHe_i = zeros(N+2*Ncp,N+2*Ncp,nT);
        R_hhe_i(:,:) = R_hh_err(1:numPaths*nT,1:numPaths*nT);
        for iT = 0:nT-1
            R_HHe_i(:,:,iT+1) = (N+2*Ncp)*F_L*R_hhe_i(iT*numPaths+1:numPaths+iT*numPaths,iT*numPaths+1:numPaths+iT*numPaths)*F_L';
        end
        
        R_X1X1 = FIcpA*blkdiag(diag(Dd1_indices),diag(diag(R_dd(1:N,1:N))))*FIcpA';
    end
    end
       %  R_PsiPsi_new1 = JH_hatIA*R_dd_tilde*JH_hatIA';
    if strcmp(est_type,'LMMSE_CP')
         y_in = reshape(y2, [N_tot nR]);
        
         if nT ==1 
         R_PsiPsi_new2 = kron(eye(nR),JF* (R_X1X1 .*  R_HHe_i(:,:) ) * JF'  );    
         elseif nT == 2 
         R_PsiPsi_new2 = kron(eye(nR),JF* (R_X1X1 .*  R_HHe_i(:,:,1) + R_X1X1 .*  R_HHe_i(:,:,2) ) * JF'  );
         elseif nT == 4
             R_psi = zeros(N_2blk,N_2blk);
             for iT = 1:nT
                R_psi = R_psi + R_X1X1 .*  R_HHe_i(:,:,iT);
             end
             R_PsiPsi_new2 = kron(eye(nR),JF* R_psi * JF'  );
         end
         R_PsiPsi_new = R_PsiPsi_new1 + R_PsiPsi_new2;
         
         N_in = N_2blk;
         XpF_in = pilots_conf.XpF_cp;
         Fbig_in = conf.Fbig_cp;
         F_L_in = conf.F_Lcp;
    elseif strcmp(est_type,'LMMSE_CPZm')
        y_in = reshape(y2, [N_tot nR]);        
        y_in = y_in(Ncp+1:end,:);        
        N_in = N;
        Fbig_in = conf.Fbig_old;
        F_L_in = conf.F_L_old;
        XpF_in = pilots_conf.XpF_old_Zm;
        R_PsiPsi_new = 0;              
    else
        y_in = reshape(y2, [N nR]);
        R_PsiPsi_new = 0;
        N_in = N;
        XpF_in = pilots_conf.XpF_old_Zm;
        Fbig_in = conf.Fbig_old;
        F_L_in = conf.F_L_old;
    end
 %   R_PsiPsi_new = calc_Rpsipsi('R_PsiPsi_CP_Dhat', nT, nR, conf, Dd_tilde, Dd1_indices);
    
    if ~strcmp(est_type,'LMMSE_CPZm')
    [estimation_out, R_hh_err] = do_Chan_Estimation_mimo(est_type, y_in, snr, N_in, nT, nR, Fbig_in, F_L_in, conf.A, XpF_in, P, numPaths, path_delays, pilots_conf.Dp, Delta_k, R_PsiPsi_new, Ncp, p.M, p.K);
    end
    
    if strcmp(est_type,'LMMSE_Zm')
        vH_lmmse_in(:,its+1) = estimation_out.vH_lmmse;
    elseif strcmp(est_type,'LMMSE_CPZm')
        [estimation_out, ~] = do_Chan_Estimation_mimo('LMMSE_Zm', y_in, snr, N_in, nT, nR, Fbig_in, F_L_in, conf.A, XpF_in, P, numPaths, path_delays, pilots_conf.Dp, Delta_k, R_PsiPsi_new, Ncp, p.M, p.K);
        
        y_in = reshape(y2, [N_tot nR]);        
        y_in = y_in(1:N,:);
        N_in = N+Ncp;
        Fbig_in = conf.Fbig_cp;
        F_L_in = conf.F_Lcp;
        XpF_in = pilots_conf.XpF_cp_2_Zm;
        R_PsiPsi_new = estimation_out;       
        
        [vH_lmmse_in(:,its+1), R_hh_err] = do_Chan_Estimation_mimo('LMMSE_CPZm', y_in, snr, N_in, nT, nR, Fbig_in, F_L_in, conf.A, XpF_in, P, numPaths, path_delays, pilots_conf.Dp, Delta_k, R_PsiPsi_new, Ncp, p.M, p.K);
                
      %  vH_lmmse_in(:,its+1) = 0.5*(vH_lmmse_in(:,its+1) + vH_lmmse);
    else
        vH_lmmse_in(:,its+1) = estimation_out;
    end    
    
    %      test = do_Chan_Estimation_mimo('LMMSE_CP', y2, snr(si), N_2blk, nT, nR, Fbig_cp, F_Lcp, 0, XpF_cp, P, numPaths, Dd, Delta_k, R_PsiPsi_new_all,Ncp);
    
    %      mean(abs(vH_lmmse_in(:,sth+1) - test).^2)
    
    %     mean(abs(vH_big - test).^2)
   
%     if (strcmp(est_type,'LMMSE_Zm') || strcmp(est_type,'LMMSE_CPZm'))
%        if (actual_snr >= 25 && its < no_PIC_iterations)
%           snr = snr+3; % (its)*snr_step+15;
%        end    
%     end 
        
    
end

wH_lmmse_in = reshape(vH_lmmse_in(:,no_PIC_iterations+1), [N nT*nR]);    
% wH_lmmse_in = reshape(conf.vH_big, [N nT*nR]);   
pilots_conf.blk = 1;
[d_hat1, ~, equalizer1, R_dd_hat1, HA1] = do_CEPIC_equalization( conf, pilots_conf, snr, wH_lmmse_in, R_hh_err, Y_all_v(:,1));

[d_hat1_soft, Heq1] = calc_soft_equichan_LMMSE(d_hat1, pilots_conf, equalizer1, R_dd_hat1, HA1, nT);

pilots_conf.blk = 2;
[d_hat2, ~, equalizer2, R_dd_hat2, HA2] = do_CEPIC_equalization( conf, pilots_conf, snr, wH_lmmse_in, R_hh_err, Y_all_v(:,2));

[d_hat2_soft, Heq2] = calc_soft_equichan_LMMSE(d_hat2, pilots_conf, equalizer2, R_dd_hat2, HA2, nT);

%d_hat = [d_hat1; d_hat2];
clear d_hat
d_hat{1} = d_hat1_soft;
d_hat{2} = d_hat2_soft;

Heq{1} = Heq1;
Heq{2} = Heq2;

vH_lmmse_out = vH_lmmse_in;


end

