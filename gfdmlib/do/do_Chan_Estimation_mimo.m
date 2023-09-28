function [ Hhat, error_cov ] = do_Chan_Estimation_mimo(method, Y, snr, N, nT, nR, Fbig, F_L, A, XpF, P, numPaths, path_delays, Dp, Delta_k, R_PsiPsi, Ncp, M, Kp, on_bins1)
%
% method:
%       1: Estimate the channel using only one pilot subsymbol
%       2: Estimate the channel through Harmonic mean of the pilot subsymbols
%       
%
%
% X:    FFT of transmit signal
% Y:    FFT of received signal
% Dp:   K*M Matrix of Pilot signals with zero on data positions
% numPaths: Number of paths in FSC situation
%
% Author: Shahab Ehsanfar
error_cov = 0;
switch method
    case 'LS'
        h_ls = (XpF' * XpF)\XpF'*Y/sqrt(N);
        
        % Vectorize the estimated channel
        VecMet = 'default';        
        %VecMet = 'Extract Channels';        
        if strcmp(VecMet,'Extract Channels')
            % tic
            h_ls_i = zeros(numPaths,nT,nR); H_ls_i = zeros(N/Delta_k,nT,nR);
            %H_i = zeros(p.M*p.K,nT,nR);
            vH_ls = []; %vH = [];
            for iR = 1:nR
                for iT = 1:nT
                    h_ls_i(:,iT,iR) = h_ls(numPaths*(iT-1)+1:numPaths*iT,iR);
                    H_ls_i(:,iT,iR) = sqrt(N)*F_L*h_ls_i(:,iT,iR);
                    %H_i(:,iT,iR) = fft(h_i(:,iT,iR),p.M*p.K);
                    
                    vH_ls = [vH_ls; H_ls_i(:,iT,iR)];
                    %vH = [vH; H_i(:,iT,iR)];
                end
            end
            % toc
        else
            % tic;
            H_ls_mat = sqrt(N)*Fbig*h_ls;
            vH_ls = H_ls_mat(:); %reshape(sqrt(N)*Fbig*h_ls,[N/Delta_k*nT*nR 1]);
            % toc
        end
        
        
        Hhat = vH_ls;
 %   Hhat = h_ls;
    
    case 'LMMSE'            
        
        R_hh = diag(P(:));
        
        VarS = 1; % Signal power
        VarN = VarS/(10^(snr/10));
        
        xp_tild = kron(eye(nR),XpF);
        
%         Sigma_yy = N*xp_tild*R_hh*xp_tild' + R_PsiPsi + VarN*eye(N/Delta_k*nR);
        Sigma_hy = R_hh*xp_tild(on_bins1,:)'.*sqrt(N);
         
%         MMSE = Sigma_hy/Sigma_yy;
       with_interference = true;
       
       if with_interference
           R_ee = R_PsiPsi(on_bins1,on_bins1)+ VarN*eye(length(on_bins1)); %eye(N*nR);
           xp_tild_spy = sparse(xp_tild(on_bins1,:));
           inv_term = xp_tild_spy'/sparse(R_ee);
       else
           xp_tild_spy = sparse(xp_tild);
           inv_term = xp_tild_spy'/VarN;
       end
        
        MMSE = (diag(1./(N*diag(R_hh))) + inv_term*xp_tild_spy )\inv_term./sqrt(N);
        Y_vectorized = Y(:);
        h_lmmse = MMSE*Y_vectorized(on_bins1);    
        error_cov = R_hh - MMSE*Sigma_hy';
       
%         R_PsiPsi(abs(R_PsiPsi)<1e-15) = 0;
%         indx1 = find(diag(R_PsiPsi),1);        
%         blk1 = R_PsiPsi(indx1:indx1+M-1,indx1:indx1+M-1) + VarN*eye(M);
%     %   blk2 = R_PsiPsi(indx1+M:indx1+2*M-1,indx1+M:indx1+2*M-1) + VarN*eye(M);        
%     %    assert( mean(abs(blk1(:) - blk2(:) ).^2) < 1e-15 ,'The two blocks are not equal for matrix inversion in LMMSE');        
%         blk1_inv = inv(blk1);        
% %         blks_str = '(1/VarN)*eye(indx1-1)';
% %         for iBlk = 1:Kp
% %             blks_str = [blks_str ',blk1_inv'];
% %         end
% %         Ree_inv = eval(['blkdiag(' blks_str ',(1/VarN)*eye(indx1-2))']);
%         %Ree_inv = blkdiag((1/VarN)*eye(indx1-1),kron(eye(Kp),blk1_inv),(1/VarN)*eye(indx1-2));
%         subRee_inv = kron(eye(Kp),blk1_inv);
%         subxp_tilde = xp_tild(indx1:end-indx1+2,:);
%         MMSE2 = (inv(N*R_hh) + (subxp_tilde'*subRee_inv)*subxp_tilde)\(subxp_tilde'*subRee_inv)./sqrt(N);
%         h_lmmse = MMSE2*Y(indx1:end-indx1+2);                            
%         
%               
%         error_cov = 0; %R_hh - MMSE2*Sigma_hy';
                  
        
  %      mean( abs( h_lmmse(:) - h_lmmse2(:) ).^2 )
        
        
        
%         H_lmmse_i = zeros(N/Delta_k,nT,nR); vH_lmmse = [];
%         for iR = 1:nR
%             for iT = 1:nT
%                 h_lmmse_i = h_lmmse( (iR-1)*nT*numPaths + (iT-1)*numPaths +(1:numPaths));
%                 H_lmmse_i(:,iT,iR) = sqrt(N)*F_L*h_lmmse_i;
%                 vH_lmmse = [vH_lmmse; H_lmmse_i(:,iT,iR)];
%             end
%         end
       % Hhat = vH_lmmse;  
        Hhat = h_lmmse;
        
    case 'HM'
        % In MIMO, HM is worse than LS
        
        Xp1F = XpF.part1;
        Xp2F = XpF.part2;
        
        h_ls_1 = (Xp1F' * Xp1F)\Xp1F'*Y/sqrt(N);
        h_ls_2 = (Xp2F' * Xp2F)\Xp2F'*Y/sqrt(N);
        
        h_ls = 0.5*(h_ls_1 + h_ls_2);
               
        % Vectorize the estimated channel
        VecMet = 'default';        
        %VecMet = 'Exctract Channels';        
        if strcmp(VecMet,'Exctract Channels')
            % tic
            h_ls_i = zeros(numPaths,nT,nR); H_ls_i = zeros(N/Delta_k,nT,nR);
            %H_i = zeros(p.M*p.K,nT,nR);
            vH_ls = []; %vH = [];
            for iR = 1:nR
                for iT = 1:nT
                    h_ls_i(:,iT,iR) = h_ls(numPaths*(iT-1)+1:numPaths*iT,iR);
                    H_ls_i(:,iT,iR) = sqrt(N)*F_L*h_ls_i(:,iT,iR);
                    %H_i(:,iT,iR) = fft(h_i(:,iT,iR),p.M*p.K);
                    
                    vH_ls = [vH_ls; H_ls_i(:,iT,iR)];
                    %vH = [vH; H_i(:,iT,iR)];
                end
            end
            % toc
        else
            % tic;
            H_ls_mat = sqrt(N)*Fbig*h_ls;
            vH_ls = H_ls_mat(:); %reshape(sqrt(N)*Fbig*h_ls,[N/Delta_k*nT*nR 1]);
            % toc
        end       
        
        Hhat = vH_ls;
        
        
    case 'LMMSE_CP'
        
        R_hh = diag(P(:));
        
        VarS = 1; % Signal power
        VarN = VarS/(10^(snr/10));
        
        xp_tild = kron(eye(nR),XpF);
        
%        Sigma_yy = N*xp_tild*R_hh*xp_tild' + R_PsiPsi + VarN*eye((N-Ncp)*nR);
        Sigma_hy = R_hh*xp_tild'.*sqrt(N);
        
%        MMSE = Sigma_hy/Sigma_yy;
        use_alternative = false;
        
   %     R_PsiPsi(abs(R_PsiPsi)<1e-1) = 0;
        
        if use_alternative
            R_ee = R_PsiPsi+ VarN*eye((N-Ncp)*nR);
            slow_part = (xp_tild'/R_ee);
            MMSE_iR = (inv(N*R_hh) + slow_part*xp_tild)\slow_part./sqrt(N);
        else
            Sigma_yy = N*xp_tild*R_hh*xp_tild' + R_PsiPsi + VarN*eye((N-Ncp)*nR);
            
            LnT = size(P,1)*size(P,2);
            indices = 1:(N-Ncp);
            MMSE_iR = Sigma_hy(1:LnT,indices) / Sigma_yy(indices,indices);            
            % MMSE = kron(eye(nR),MMSE_iR);            
            
            % MMSE = Sigma_hy/Sigma_yy;
        end
        % h_lmmse = MMSE*Y(:);
        h_lmmse_w = zeros(LnT,nR);
        for iR = 1:nR
            h_lmmse_w(:,iR) = MMSE_iR*Y(:,iR);
        end
        h_lmmse = h_lmmse_w(:);
        
        R_hh_hat = kron(eye(nR),MMSE_iR*Sigma_hy(1:LnT,indices)'); % MMSE*Sigma_hy';
        error_cov = R_hh - R_hh_hat;
        
        H_lmmse_i = zeros(N-2*Ncp,nT,nR); vH_lmmse = [];
        h_lmmse_i = zeros(path_delays(end),1);
        for iR = 1:nR
            for iT = 1:nT
                h_lmmse_i(path_delays) = h_lmmse( (iR-1)*nT*numPaths + (iT-1)*numPaths +(1:numPaths));
                H_lmmse_i(:,iT,iR) = fft(h_lmmse_i,N-2*Ncp);%sqrt(N-Ncp)*F_L*h_lmmse_i;
                vH_lmmse = [vH_lmmse; H_lmmse_i(:,iT,iR)];
            end
        end
        Hhat = vH_lmmse;  
    %    Hhat = fft(h_lmmse,N/2);
    
    case 'LMMSE_CPZ'
        
        R_hh = diag(P(:));
        
        VarS = 1; % Signal power
        VarN = VarS/(10^(snr/10));
        
        xp_tild = kron(eye(nR),XpF);
        
%        Sigma_yy = N*xp_tild*R_hh*xp_tild' + R_PsiPsi + VarN*eye((N-Ncp)*nR);
        Sigma_hy = R_hh*xp_tild'.*sqrt(N+Ncp);
        
%        MMSE = Sigma_hy/Sigma_yy;
        use_alternative = false;
        
        if use_alternative
            R_ee = R_PsiPsi+ VarN*eye((N-Ncp)*nR);
            slow_part = (xp_tild'/R_ee);
            MMSE2 = (inv(N*R_hh) + slow_part*xp_tild)\slow_part./sqrt(N+Ncp);
        else
            Sigma_yy = (N+Ncp)*xp_tild*R_hh*xp_tild' + R_PsiPsi + VarN*eye((N-Ncp)*nR);
            
            LnT = size(P,1)*size(P,2);
            indices = 1:(N-Ncp);
            MMSE_iR = Sigma_hy(1:LnT,indices) / Sigma_yy(indices,indices);            
            % MMSE = kron(eye(nR),MMSE_iR);            
            
            % MMSE = Sigma_hy/Sigma_yy;
        end
        % h_lmmse = MMSE*Y(:);
        h_lmmse_w = zeros(LnT,nR);
        for iR = 1:nR
            h_lmmse_w(:,iR) = MMSE_iR*Y(:,iR);
        end
        h_lmmse = h_lmmse_w(:);
        
        R_hh_hat = kron(eye(nR),MMSE_iR*Sigma_hy(1:LnT,indices)'); % MMSE*Sigma_hy';
        error_cov = R_hh - R_hh_hat;
        
        H_lmmse_i = zeros(N-Ncp,nT,nR); vH_lmmse = [];
        h_lmmse_i = zeros(path_delays(end),1);
        for iR = 1:nR
            for iT = 1:nT
                h_lmmse_i(path_delays) = h_lmmse( (iR-1)*nT*numPaths + (iT-1)*numPaths +(1:numPaths));
                H_lmmse_i(:,iT,iR) = fft(h_lmmse_i,N-Ncp);%sqrt(N-Ncp)*F_L*h_lmmse_i;
                vH_lmmse = [vH_lmmse; H_lmmse_i(:,iT,iR)];
            end
        end
        Hhat = vH_lmmse;  
    %    Hhat = fft(h_lmmse,N/2);
    
    case 'LowComp_LMMSE_Z'
        
        R_hh = diag(P(:));
        
        VarS = 1; % Signal power
        VarN = VarS/(10^(snr/10));
        
        xp_tild = kron(eye(nR),XpF);
        
        Sigma_yy = N*xp_tild*R_hh*xp_tild' + R_PsiPsi + VarN*eye(N/Delta_k*nR);
        Sigma_hy = R_hh*xp_tild'.*sqrt(N);
         
        MMSE = Sigma_hy/Sigma_yy;
       
%        R_ee = R_PsiPsi+ VarN*eye(N/Delta_k*nR);
%        MMSE = (inv(N*R_hh) + (xp_tild'/R_ee)*xp_tild)\(xp_tild'/R_ee)./sqrt(N);
        h_lmmse = MMSE*Y(:);    
        error_cov = R_hh - MMSE*Sigma_hy';
    
    case 'LMMSE_CPZm' % Ignores the interference
        F_M = dftmtx(M)/sqrt(M);
        h_lmmse_m = zeros(numPaths*nT*nR,M);
        error_cov_m = zeros(numPaths*nT,numPaths*nT,M);
        variances_m_iT = zeros(numPaths*nT,M);
%         m = 1;
%             Z_m = kron(F_M(m,:),eye(Kp));
%             Zy_cp = Z_m*Y;
%             
%             P_iR1 = P(:,:,1);
%             R_hh = diag(P_iR1(:));
%             
             VarS = 1; % Signal power
             VarN = VarS/(10^(snr/10));
%             
%             xp_tild_iR = sqrt(N+Ncp)*XpF(:,:,m);
%             %xp_tild = kron(eye(nR),XpF(:,:,m));
%             %xp_tild = kron(eye(nR),([XpF(:,:,1); XpF(:,:,2);XpF(:,:,3);XpF(:,:,4);XpF(:,:,5);XpF(:,:,6);XpF(:,:,7)]));
%             
%             %        Sigma_yy = N*xp_tild*R_hh*xp_tild' + R_PsiPsi + VarN*eye((N-Ncp)*nR);
%             %Sigma_hy = R_hh*xp_tild'.*sqrt(N+Ncp);
%             Sigma_hy = R_hh*xp_tild_iR';
%             
%             %        MMSE = Sigma_hy/Sigma_yy;
%             use_alternative = true;
%             
%             if use_alternative
%                 R_ee =   VarN*eye((Kp)); % R_PsiPsi(:,:,m) +
%                 slow_part = (xp_tild_iR'/R_ee);
%                 MMSE_iR = (inv(R_hh) + slow_part*xp_tild_iR)\slow_part;
%           %     diagonaldominant = (inv((N+Ncp)*R_hh) + slow_part*xp_tild);
%           %     MMSE_iR = (diag(1./diag(diagonaldominant)))*slow_part./sqrt(N+Ncp);
%                 LnT = size(P,1)*size(P,2);
%                 indices = 1:(Kp);
%             else
%                 % IGNORE THE INTERFERENCE TERM
%                 Sigma_yy = (N+Ncp)*xp_tild*R_hh*xp_tild' + VarN*eye((Kp)*nR);
%                 
%                 LnT = size(P,1)*size(P,2);
%                 indices = 1:(Kp);
%                 MMSE_iR = Sigma_hy(1:LnT,indices) / Sigma_yy(indices,indices);
%                 % MMSE = kron(eye(nR),MMSE_iR);
%                 
%                 % MMSE = Sigma_hy/Sigma_yy;
%             end
%             % h_lmmse = MMSE*Y(:);
%             h_lmmse_w = zeros(LnT,nR);
%             for iR = 1:nR
%                 h_lmmse_w(:,iR) = MMSE_iR*Zy_cp(:,iR);
%             end
%             h_lmmse_m(:,m) = h_lmmse_w(:);
%             error_cov_m(:,:,m) = inv((inv(R_hh) + slow_part*xp_tild_iR/VarN ));          
%           %  R_hh_hat = kron(eye(nR),MMSE_iR*Sigma_hy(1:LnT,indices)'); % MMSE*Sigma_hy';
%           %  error_cov_m(:,:,m) = kron(eye(nR),R_hh) - R_hh_hat;
%             
%         %    mean(abs( h_lmmse_w(:) - A(:) ).^2)
            PriorEst = R_PsiPsi;
            m = 1;
            h_lmmse_w = PriorEst.h_lmmse;
            error_cov_m(:,:,m) = PriorEst.error_cov;
            
    %       mean(abs( h_lmmse_w(:) - A(:) ).^2)
            
            while (m< M+1)         
            
            if m > 0
            Z_m = kron(F_M(m,:),eye(Kp));
            Zy_cp = Z_m*Y;
            xp_tild_iR = sqrt(N+Ncp)*XpF(:,:,m);
            
            R_y_tilde = xp_tild_iR*error_cov_m(:,:,m)*xp_tild_iR' + VarN*eye(Kp);
            R_hy_tilde = error_cov_m(:,:,m)*xp_tild_iR';
            
            estimator = R_hy_tilde/R_y_tilde;
            h_new = h_lmmse_w + estimator*(Zy_cp - xp_tild_iR*h_lmmse_w);
            error_cov_m(:,:,m+1) = (eye(numPaths*nT) - estimator*xp_tild_iR)*error_cov_m(:,:,m);
            else
                error_cov_m(:,:,m+1) = error_cov_m(:,:,m);
                h_new = h_lmmse_w;
            end
            
            h_lmmse_m(:,m+1) = h_new(:);
            
%            Er = mean(abs( h_new(:) - A(:) ).^2);
%            [m Er]
            
            
            m = m+1;
            h_lmmse_w = h_new;
            end          
                
         error_cov = kron(eye(nR),error_cov_m(:,:,end));
         h_lmmse = h_lmmse_m(:,end);
        
        
         H_lmmse_i = zeros(N-Ncp,nT,nR); vH_lmmse = [];
         h_lmmse_i = zeros(path_delays(end),1);
         for iR = 1:nR
             for iT = 1:nT
                 h_lmmse_i(path_delays) = h_lmmse( (iR-1)*nT*numPaths + (iT-1)*numPaths +(1:numPaths));
                 H_lmmse_i(:,iT,iR) = fft(h_lmmse_i,N-Ncp);%sqrt(N-Ncp)*F_L*h_lmmse_i;
                 vH_lmmse = [vH_lmmse; H_lmmse_i(:,iT,iR)];
             end
         end
         Hhat = vH_lmmse;
         %    Hhat = fft(h_lmmse,N/2);
         
    case 'LMMSE_Zm' % Ignores the interference
        F_M = dftmtx(M)/sqrt(M);
      %  h_lmmse= zeros(numPaths*nT*nR,M);
      %  for m_outer = 1:M
        h_lmmse_m = zeros(numPaths*nT*nR,M);
        error_cov_m = zeros(numPaths*nT,numPaths*nT,M);
        variances_m_iT = zeros(numPaths*nT,M);
      %  Z = kron(F_M,eye(Kp));
        m = 1;
            indx = m; %mod(m+m_outer-1,7)+1;
            Z_m = kron(F_M(m,:),eye(Kp));
            Zy_cp = Z_m*Y;
            
            P_iR1 = P(:,:,1);
            R_hh = diag(P_iR1(:));
            
            VarS = 1; % Signal power
            VarN = VarS/(10^(snr/10));
            
            xp_tild_iR = sqrt(N)*Z_m*XpF(:,:,indx);
       %     xp_tild_iR_m = Z_m*XpF(:,:,indx);
            %xp_tild = kron(eye(nR),XpF(:,:,m));
            %xp_tild = kron(eye(nR),([XpF(:,:,1); XpF(:,:,2);XpF(:,:,3);XpF(:,:,4);XpF(:,:,5);XpF(:,:,6);XpF(:,:,7)]));
            
            %        Sigma_yy = N*xp_tild*R_hh*xp_tild' + R_PsiPsi + VarN*eye((N-Ncp)*nR);
            Sigma_hy = R_hh*xp_tild_iR'.*sqrt(N);
            
            %        MMSE = Sigma_hy/Sigma_yy;
            use_alternative = true;
            
            if use_alternative
                % 
               slow_part = (xp_tild_iR');
               MMSE_iR = (inv(R_hh) + slow_part*xp_tild_iR/VarN)\xp_tild_iR'/(VarN);%slow_part*xp_tild_iR
%                 R_ee = R_PsiPsi(:,:,m) + VarN*eye((Kp)); % 
%                 slow_part = (xp_tild_iR')/R_ee;
%                 MMSE_iR = (inv(R_hh) + slow_part*xp_tild_iR)\slow_part;
               
                LnT = size(P,1)*size(P,2);
                indices = 1:(Kp);
            else
                % IGNORE THE INTERFERENCE TERM
                Sigma_yy = (N)*xp_tild_iR*R_hh*xp_tild_iR' + VarN*eye((N)*nR);
               %Sigma_yy = (N)*xp_tild_iR*R_hh*xp_tild_iR' + R_PsiPsi(:,:,m) + VarN*eye((N)*nR);
                
                LnT = size(P,1)*size(P,2);
                indices = 1:(Kp);
                MMSE_iR = Sigma_hy(1:LnT,indices) / Sigma_yy(indices,indices);
                % MMSE = kron(eye(nR),MMSE_iR);
                
                % MMSE = Sigma_hy/Sigma_yy;
            end
            % h_lmmse = MMSE*Y(:);
            h_lmmse_w = zeros(LnT,nR);
            for iR = 1:nR
                h_lmmse_w(:,iR) = MMSE_iR*Zy_cp(:,iR);
            end
            h_lmmse_m(:,m) = h_lmmse_w(:);
            
            error_cov_m(:,:,m) = inv((inv(R_hh) + slow_part*xp_tild_iR/VarN ));
%             error_cov_m(:,:,m) = inv((inv(R_hh) + slow_part*xp_tild_iR ));
           
            %  R_hh_hat = kron(eye(nR),MMSE_iR*Sigma_hy(1:LnT,indices)'); % MMSE*Sigma_hy';
          %  error_cov_m(:,:,m) = kron(eye(nR),R_hh) - R_hh_hat;
          %  error_cov_m(:,:,m) = R_hh - MMSE_iR*Sigma_hy(1:LnT,indices)';
       %     mean(abs( h_lmmse_w(:) - A(:) ).^2)
            while (m< M)         
            
            Z_m = kron(F_M(m+1,:),eye(Kp));
            Zy_cp = Z_m*Y;
            xp_tild_iR = sqrt(N)*Z_m*XpF(:,:,m+1);
            
            R_y_tilde = xp_tild_iR*error_cov_m(:,:,m)*xp_tild_iR' + VarN*eye(Kp);
%            R_y_tilde = xp_tild_iR*error_cov_m(:,:,m)*xp_tild_iR' + R_PsiPsi(:,:,m) + VarN*eye(Kp);            
            R_hy_tilde = error_cov_m(:,:,m)*xp_tild_iR';
           
            
            estimator = R_hy_tilde/R_y_tilde; 
            h_new = h_lmmse_w + estimator*(Zy_cp - xp_tild_iR*h_lmmse_w);
            h_lmmse_m(:,m+1) = h_new(:);
            
    %        mean(abs( h_new(:) - A(:) ).^2)
            
            error_cov_m(:,:,m+1) = (eye(numPaths*nT) - estimator*xp_tild_iR)*error_cov_m(:,:,m);
            m = m+1;
            h_lmmse_w = h_new;
            end
            
        error_cov = kron(eye(nR),error_cov_m(:,:,end));    
        
        
        h_lmmse = h_lmmse_m(:,end);%(:,m_outer)
    %    end
        
         H_lmmse_i = zeros(N,nT,nR); vH_lmmse = [];
         h_lmmse_i = zeros(path_delays(end),1);
         for iR = 1:nR
             for iT = 1:nT
                 h_lmmse_i(path_delays) = h_lmmse( (iR-1)*nT*numPaths + (iT-1)*numPaths +(1:numPaths));
                 H_lmmse_i(:,iT,iR) = fft(h_lmmse_i,N);%sqrt(N-Ncp)*F_L*h_lmmse_i;
                 vH_lmmse = [vH_lmmse; H_lmmse_i(:,iT,iR)];
             end
         end
         Hhat.vH_lmmse = vH_lmmse;
         Hhat.h_lmmse = h_lmmse_w;
         Hhat.error_cov = error_cov_m(:,:,end);
         %    Hhat = fft(h_lmmse,N/2);    
        
    case 'LS_Zm'
        
        
        F_M = dftmtx(M)/sqrt(M);
        if nT == 1
            h_test = zeros(numPaths,M);
        else
            h_test = zeros(nR*numPaths,nT,M);
        end
        norms = zeros(1,M);
        for m = 1:M            
            
            Z_m = kron(F_M(m,:),eye(Kp));
            
            Ym = Z_m*Y;
            
            for iT = 1:nT
            dp = Dp(:,:,iT);
            dp = dp(:);
            Zdp = Z_m*dp;
            Zdd = 1;%trace(Z_m*A*diag([zeros(p.K,1); ones(p.K*p.Md,1); zeros(p.K,1)])*A*Z_m');
                        
            norms(m) = (Zdp'*Zdp)/abs(Zdd);
                        
%             A1 = A(:,1:Kp);
%             
%             xp_z = A1*Z_m*dp;
            end
%             Qp = Z_m*F'*diag(F*xp_z)*F_L;
%             norm2(m) = trace(Qp'*Qp);
            
%             Xp_new = Z_m*F'*diag(F*xp_z)*F_L;
            
            if nT == 1
                h_test(:,m) = norms(m)*((XpF(:,:,m)' * XpF(:,:,m))\XpF(:,:,m)'*Ym/Kp);%(Xp_new'*Xp_new)\Xp_new'*Ym/Kp;%(XpF' * XpF)\XpF'*Y/sqrt(N);
            else
                h_test(:,:,m) = norms(m)*((XpF(:,:,m)' * XpF(:,:,m))\XpF(:,:,m)'*Ym/Kp);
            end
            
        end        
        
        weights = repmat(norms,[numPaths 1]);
        
        if nT == 1
            h_ls = sum(h_test,2)./sum(weights,2);
        else
            h_ls = sum(h_test,3)./sum(weights(1,:),2);
        end

        % Vectorize the estimated channel
        VecMet = 'default';        
        %VecMet = 'Extract Channels';        
        if strcmp(VecMet,'Extract Channels')
            % tic
            h_ls_i = zeros(numPaths,nT,nR); H_ls_i = zeros(N/Delta_k,nT,nR);
            %H_i = zeros(p.M*p.K,nT,nR);
            vH_ls = []; %vH = [];
            for iR = 1:nR
                for iT = 1:nT
                    h_ls_i(:,iT,iR) = h_ls(numPaths*(iT-1)+1:numPaths*iT,iR);
                    H_ls_i(:,iT,iR) = sqrt(N)*F_L*h_ls_i(:,iT,iR);
                    %H_i(:,iT,iR) = fft(h_i(:,iT,iR),p.M*p.K);
                    
                    vH_ls = [vH_ls; H_ls_i(:,iT,iR)];
                    %vH = [vH; H_i(:,iT,iR)];
                end
            end
            % toc
        else
            % tic;
            H_ls_mat = sqrt(N)*Fbig*h_ls;
            vH_ls = H_ls_mat(:); %reshape(sqrt(N)*Fbig*h_ls,[N/Delta_k*nT*nR 1]);
            % toc
        end
        
        
        Hhat = vH_ls;
    
    case 'LowComp_LMMSE_CP_DRAFT'
        
        
        F_M = dftmtx(p.M)/sqrt(p.M);
        Z = kron(F_M,eye(p.K));
        
        h_circ = F_cpblk'*diag(fft(h_gains,N+2*Ncp))*F_cpblk;
        %J = [zeros(N+Ncp,Ncp) eye(N+Ncp)];
        %J = [zeros(N,2*Ncp) eye(N)];
        J_cp = [zeros(N,Ncp) eye(N) zeros(N,Ncp)];
        I_cp;
        
        JHIA = J_cp*h_circ*IcpA;%*[zeros(N);eye(N)];
        %JHIA_cp = JHIA(1:N,:);
        %JHIA_normal = JHIA(Ncp+1:end,:);
        Z_2blk = kron(eye(2),Z);
        
        Zycp = Z*ycp(1:N);
        xp_z = F_cpblk*IcpA*[zeros(N,1); Dp(:)];
        ZXp = Z*J_cp*F_cpblk'*diag(xp_z)*Fbig_Lcpblk;
        
        h_hat_cp_big = (ZXp'*ZXp)\ZXp'*Zycp/sqrt(N+2*Ncp);
        
        
        
        
        
               
        h_test = zeros(numPaths,p.M);
        for m = 1:p.M
            if m == 44 || m == 54 %|| m == 3 || m == 6
                
            else
                Z_m = kron(F_M(m,:),eye(p.K));
                
                Ym = Z_m*ycp(1:N);
                
                %         dp = [zeros(N,1); Dp(:)];
                %         Zdp = kron(eye(2),Z_m)*dp;
                %
                %         dp = Dp(:);
                %         Zdp = Z_m*dp;
                %         figure; plot(abs(Zdp))
                %         ZJHIAZ = Z*JHIA*Z_2blk';
                %         ZJHIAZm = ZJHIAZ(1:p.K,337:336+p.K);
                %         ytest = ZJHIAZm*Zdp;
                %
                %         ymtest = Z_m*J_cp*F_cpblk'*diag(fft(h_gains,N+2*Ncp))*F_cpblk*IcpA*[zeros(N,1); Dp(:)]; % *kron(eye(2),Z_m')*Zdp;
                xp_z = F_cpblk*IcpA*[zeros(N,1); Dp(:)];
                %        ymtest = Z_m*J_cp*F_cpblk'*diag(xp_z)*Fbig_Lcpblk*h_gains*sqrt(432);
                
                
                m
                %norm = Zdp'*Zdp
                
                %Icp_A1 = IcpA(:,0*p.K+1:2*p.K);
                
                %xp_z = Icp_A1*Zdp;
                %xp_z_cp = [zeros(Ncp,1); xp_z(end-Ncp+1:end); xp_z];
                
                Xp_new = Z_m*J_cp*F_cpblk'*diag(xp_z)*Fbig_Lcpblk;
                %Xp_new = Z_m*F'*diag(F*xp_z)*F_L;
                
                h_test(:,m) = (Xp_new'*Xp_new)\Xp_new'*Ym/sqrt(432);
            end
        end
        %figure; plot(abs(h_test)) %%%% WORKS WORKS WORKS
        
        h_test_all = mean(h_test,2);
        
        figure(102);
        hold on; plot(abs(h_test_all),'b-+')
        
        hold on; plot(abs(h_gains),'k')
        
        mean(abs(h_test_all - h_gains).^2)
        
        
        figure; plot(abs(h_hat_cp_big),'b-+')
        hold on; plot(abs(h_gains),'k')
 
        mean(abs(h_hat_cp_big - h_gains).^2)
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
  
    
end

