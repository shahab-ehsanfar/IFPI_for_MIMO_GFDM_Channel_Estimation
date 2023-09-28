function [ H_eq, Hhat ] = do_Chan_Estimation_mimo_ser(method, Y, snr, N, nT, nR, Fbig, F_L, FA, XpF, P, numPaths, Dd, Delta_k, R_PsiPsi)
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
        h_ls = (XpF' * XpF)\XpF'*Y/sqrt(N);
                         
       % H_ls = sqrt(N)*Fbig*h_ls;
        
        H_eq = zeros(N*nT,N*nR);
        H_i = zeros(N,nT,nR);
        for iR = 1:nR            
            for iT = 1:nT
                H_i(:,iT,iR) = fft(h_ls((iT-1)*numPaths+1:(iT-1)*numPaths+numPaths,iR),N);
                for n = 1:N
                    H_eq((iR-1)*N+n,(iT-1)*N+n) = H_i(n,iT,iR);
                end
            end
        end
        Hhat = H_i;
    
    case 'LMMSE'            
        
        R_hh = diag(P(:).^2);
        
        VarS = 1; % Signal power
        VarN = VarS/(10^(snr/10));
        
        xp_tild = kron(eye(nR),XpF);
        
        Sigma_yy = N*xp_tild*R_hh*xp_tild' + R_PsiPsi + VarN*eye(N/Delta_k*nR);
        Sigma_hy = R_hh*xp_tild'.*sqrt(N);
        
        h_lmmse = Sigma_hy/Sigma_yy*Y(:);
        
        H_lmmse_i = zeros(N,nT,nR); %vH_lmmse = [];
        H_eq = zeros(N*nT,N*nR);
        for iR = 1:nR
            for iT = 1:nT
                h_lmmse_i = h_lmmse( (iR-1)*nT*numPaths + (iT-1)*numPaths +(1:numPaths));
                %H_lmmse_i(:,iT,iR) = sqrt(N)*F_L*h_lmmse_i;
                H_lmmse_i(:,iT,iR) = fft(h_lmmse_i,N);
                %vH_lmmse = [vH_lmmse; H_lmmse_i(:,iT,iR)];
                for n = 1:N
                    H_eq((iR-1)*N+n,(iT-1)*N+n) = H_lmmse_i(n,iT,iR);
                end
                %H_eq((iR-1)*N+(1:N),(iT-1)*N+(1:N)) = diag(H_lmmse_i(:,iT,iR));
            end
        end
        %H_lmmse = reshape(H_lmmse_i, [N*nT nR]);
        
        Hhat = H_lmmse_i;        
    
end

