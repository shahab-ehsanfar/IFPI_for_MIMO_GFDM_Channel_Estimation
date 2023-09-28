function [R_PsiPsi, aEi, aEn ] = calc_ce_mse_mimo(method, Y, snr, N, nT, nR, Fbig, F_L, FA, XpF, P, numPaths, Dd, Delta_k, R_PsiPsi)
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

switch method
    case 'LS'
        Q_ls = (XpF' * XpF)\XpF'; %*Y/sqrt(N);
        
        I_kron_QhQ = (1/(Delta_k))*kron(eye(nR),Q_ls'*Q_ls);
        
        aEi = real(trace(I_kron_QhQ*R_PsiPsi)/(N/Delta_k*nR*nT));
        
        VarS = 1; % Signal power
        VarN = VarS/(10^(snr/10));
        
        aEn = real(VarN*trace(I_kron_QhQ)/(N/Delta_k*nR*nT));
        
        R_PsiPsi = [];
    
    case 'LMMSE'                   
        R_hh = diag(P(:));
        
        VarS = 1; % Signal power
        VarN = VarS/(10^(snr/10));
        
        xp_tild = kron(eye(nR),XpF);
        
        Sigma_yy = N*xp_tild*R_hh*xp_tild' + R_PsiPsi + VarN*eye(N/Delta_k*nR);
        Sigma_hy = R_hh*xp_tild'.*sqrt(N);
        
        I_kron_F = kron(eye(nR),Fbig);
        
        R_HH = N*I_kron_F*R_hh*I_kron_F';
        
%        % EXACT COMPUTATION        
%  %      tic
%        aEn = trace(R_HH - N*I_kron_F*Sigma_hy/Sigma_yy*Sigma_hy'*I_kron_F')/(N/Delta_k*nR*nT);
%  %      toc
        
%       % EFFICIENT COMPUTATION
%        tic
        sz = size(Sigma_yy,2)/nR;
        sz2 = size(Sigma_hy,1)/nR;
        Q = zeros(sz*nR*nT);
        for iR = 0:nR-1
           Q(1+iR*sz*nT:sz*nT+iR*sz*nT,1+iR*sz*nT:sz*nT+iR*sz*nT) = Fbig*(Sigma_hy(1+iR*sz2:sz2+iR*sz2,1+iR*sz:sz+iR*sz)/Sigma_yy(1+iR*sz:sz+iR*sz,1+iR*sz:sz+iR*sz))*Sigma_hy(1+iR*sz2:sz2+iR*sz2,1+iR*sz:sz+iR*sz)'*Fbig';
        end        
        aEn = trace(R_HH - N*Q)/(N/Delta_k*nR*nT);
%        toc
%        assert( sum(sum(Q - I_kron_F*Sigma_hy/Sigma_yy*Sigma_hy'*I_kron_F')) < 1e-8)
        
        aEi = []; R_PsiPsi = [];
        
    case 'R_PsiPsi'
        % initializations
        sigma2d = zeros(size(FA,2),nT);
        R_XdXd = zeros(N,N,nT);
        R_PsiPsi_i = zeros(N,N,nT,nR);
        R_hh_i = zeros(numPaths,numPaths,nT,nR);
        Upsilon = zeros(N,N,nT,nR);
        R_PsiPsi_iR = zeros(N,N,nR);
        
%         sigma2d = zeros(size(FA,2),nT);
%         R_XdXd = zeros(N/Delta_k,N/Delta_k,nT);
%         R_PsiPsi_i = zeros(N/Delta_k,N/Delta_k,nT,nR);
%         R_hh_i = zeros(numPaths,numPaths,nT,nR);
%         Upsilon = zeros(N/Delta_k,N/Delta_k,nT,nR);
%         R_PsiPsi_iR = zeros(N/Delta_k,N/Delta_k,nR);
        
        % interference calculation
        for iR = 1:nR
            for iT = 1:nT
                R_hh_i(:,:,iT,iR) = diag(P(:,iT,iR));
                Upsilon(:,:,iT,iR) = N*F_L*R_hh_i(:,:,iT,iR)*F_L';
                
                if size(Dd,2) ~= 1
                Dd_i = Dd(:,:,iT);
                sigma2d(:,iT) = (abs(Dd_i(:)) > 0) + 0;
                else
                    sigma2d(:,iT) = Dd;
                end
                
                R_XdXd(:,:,iT) = FA*diag(sigma2d(:,iT))*FA';
                
                R_PsiPsi_i(:,:,iT,iR) = Upsilon(:,:,iT,iR) .* R_XdXd(:,:,iT);
                R_PsiPsi_iR(:,:,iR) = R_PsiPsi_iR(:,:,iR) + R_PsiPsi_i(:,:,iT,iR);
            end
        end
        
        % Generate the large block diagonal R_PsiPsi from the individual
        % blocks
        R_PsiPsi_str = char(zeros(nR,17 + numel(num2str(nR)) ));
        for iR = 1:nR
            R_PsiPsi_str(iR,1:17+numel(num2str(iR))) = sprintf('R_PsiPsi_iR(:,:,%d)',iR);
            
            if iR == 1
                R_PsiPsi_blks = R_PsiPsi_str(iR,1:17+numel(num2str(iR)));
            else
                R_PsiPsi_blks = [R_PsiPsi_blks ',' R_PsiPsi_str(iR,1:17+numel(num2str(iR)))];
            end
        end
        R_PsiPsi = eval(['blkdiag(' R_PsiPsi_blks ');']);
        
        aEi = []; aEn = [];
        
    case 'LMMSEinterf'
        R_hh = diag(P(:).^2);
        
        VarS = 1; % Signal power
        VarN = VarS/(10^(snr/10));
        
        xp_tild = kron(eye(nR),XpF);
        
        Sigma_yy = N*xp_tild*R_hh*xp_tild' + R_PsiPsi;
        Sigma_hy = R_hh*xp_tild'.*sqrt(N);
        
        I_kron_F = kron(eye(nR),Fbig);
        
        R_HH = N*I_kron_F*R_hh*I_kron_F';
        
        %        % EXACT COMPUTATION
        %  %      tic
        %        aEn = trace(R_HH - N*I_kron_F*Sigma_hy/Sigma_yy*Sigma_hy'*I_kron_F')/(N/Delta_k*nR*nT);
        %  %      toc
        
        %       % EFFICIENT COMPUTATION
        %        tic
        sz = size(Sigma_yy,2)/nR;
        sz2 = size(Sigma_hy,1)/nR;
        Q = zeros(sz*nR*nT);
        for iR = 0:nR-1
            Q(1+iR*sz*nT:sz*nT+iR*sz*nT,1+iR*sz*nT:sz*nT+iR*sz*nT) = Fbig*(Sigma_hy(1+iR*sz2:sz2+iR*sz2,1+iR*sz:sz+iR*sz)/Sigma_yy(1+iR*sz:sz+iR*sz,1+iR*sz:sz+iR*sz))*Sigma_hy(1+iR*sz2:sz2+iR*sz2,1+iR*sz:sz+iR*sz)'*Fbig';
        end
        aEn = trace(R_HH - N*Q)/(N/Delta_k*nR*nT);
        %        toc
        %        assert( sum(sum(Q - I_kron_F*Sigma_hy/Sigma_yy*Sigma_hy'*I_kron_F')) < 1e-8)
        
        aEi = []; R_PsiPsi = [];
        
    case 'LMMSEnoise'
        R_hh = diag(P(:).^2);
        
        VarS = 1; % Signal power
        VarN = VarS/(10^(snr/10));
        
        xp_tild = kron(eye(nR),XpF);
        
        Sigma_yy = N*xp_tild*R_hh*xp_tild' + VarN*eye(N/Delta_k*nR);
        Sigma_hy = R_hh*xp_tild'.*sqrt(N);
        
        I_kron_F = kron(eye(nR),Fbig);
        
        R_HH = N*I_kron_F*R_hh*I_kron_F';
        
        %        % EXACT COMPUTATION
        %  %      tic
        %        aEn = trace(R_HH - N*I_kron_F*Sigma_hy/Sigma_yy*Sigma_hy'*I_kron_F')/(N/Delta_k*nR*nT);
        %  %      toc
        
        %       % EFFICIENT COMPUTATION
        %        tic
        sz = size(Sigma_yy,2)/nR;
        sz2 = size(Sigma_hy,1)/nR;
        Q = zeros(sz*nR*nT);
        for iR = 0:nR-1
            Q(1+iR*sz*nT:sz*nT+iR*sz*nT,1+iR*sz*nT:sz*nT+iR*sz*nT) = Fbig*(Sigma_hy(1+iR*sz2:sz2+iR*sz2,1+iR*sz:sz+iR*sz)/Sigma_yy(1+iR*sz:sz+iR*sz,1+iR*sz:sz+iR*sz))*Sigma_hy(1+iR*sz2:sz2+iR*sz2,1+iR*sz:sz+iR*sz)'*Fbig';
        end
        aEn = trace(R_HH - N*Q)/(N/Delta_k*nR*nT);
        %        toc
        %        assert( sum(sum(Q - I_kron_F*Sigma_hy/Sigma_yy*Sigma_hy'*I_kron_F')) < 1e-8)
        
        aEi = []; R_PsiPsi = [];
        
        case 'R_PsiPsi_CP'
        % initializations
        Nd = size(FA,2);
        sigma2d = zeros(Nd,nT);
        R_XdXd = zeros(N,N,nT);
        R_PsiPsi_i = zeros(N,N,nT,nR);
        R_hh_i = zeros(numPaths,numPaths,nT,nR);
        Upsilon = zeros(N,N,nT,nR);
        R_PsiPsi_iR = zeros(N,N,nR);
        I_cp = eye(Nd);
        cp_length= size(Dd,1);
        I_cp = [I_cp(end-cp_length+1:end,:); I_cp];
        
        % interference calculation
        for iR = 1:nR
            for iT = 1:nT
                R_hh_i(:,:,iT,iR) = diag(P(:,iT,iR));
                Upsilon(:,:,iT,iR) = N*F_L*R_hh_i(:,:,iT,iR)*F_L';
                
                Dd_i = Dd(:,:,iT);
                sigma2d(:,iT) = (abs(Dd_i(:)) > 0) + 0;
                R_XdXd(:,:,iT) = I_cp*FA*diag(sigma2d(:,iT))*FA'*I_cp';
                
                R_PsiPsi_i(:,:,iT,iR) = Upsilon(:,:,iT,iR) .* R_XdXd(:,:,iT);
                R_PsiPsi_iR(:,:,iR) = R_PsiPsi_iR(:,:,iR) + R_PsiPsi_i(:,:,iT,iR);
            end
        end
        
        % Generate the large block diagonal R_PsiPsi from the individual
        % blocks
        R_PsiPsi_str = char(zeros(nR,17 + numel(num2str(nR)) ));
        for iR = 1:nR
            R_PsiPsi_str(iR,1:17+numel(num2str(iR))) = sprintf('R_PsiPsi_iR(:,:,%d)',iR);
            
            if iR == 1
                R_PsiPsi_blks = R_PsiPsi_str(iR,1:17+numel(num2str(iR)));
            else
                R_PsiPsi_blks = [R_PsiPsi_blks ',' R_PsiPsi_str(iR,1:17+numel(num2str(iR)))];
            end
        end
        R_PsiPsi = eval(['blkdiag(' R_PsiPsi_blks ');']);
        
        aEi = []; aEn = [];
                
end

