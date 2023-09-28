function [ Heq, Hhat ] = do_Chan_Estimation_mimo_precode(method, p, Y, numPaths, snr, N, nT, nR, Fp, Xp_iT, h, R_HH) %Fbig, F_L, FA, XpF, P, numPaths, Dd, Delta_k, R_PsiPsi)
%
%
% Author: Shahab Ehsanfar

switch method
    case 'LS'
        y = zeros(N,nR);
        Yp = zeros(p.Kon/p.Delta_k,nT,nR);
        H_ls_i = zeros(N,nT,nR);
        Heq = zeros(N*nT,N*nR);
    %    Heq2 = zeros(N*nT,N*nR);
              
        for iR = 1:nR
            % Receive signal in time
            y(:,iR) = ifft_u(Y(:,iR));
           % y(:,iR) = awgn(ifft_u(Y(:,iR)),snr);        
            for iT = 1:nT
                Yp(:,iT,iR) = Fp(:,:,iT)*y(:,iR);
                
                % individual SISO LS estimations
                H_ls_temp = Yp(:,iT,iR)./Xp_iT(:,iT);
                
                % FFT based interpolation through zero padding
                h_ls_temp = p.Delta_k*p.M*Fp(:,:,iT)'*H_ls_temp/sqrt(N);                
                h_ls_temp(numPaths+1:end) = 0; 
                H_ls_i(:,iT,iR) = fft(h_ls_temp); 
                
                % vector to matrix
                for n = 1:N
                    Heq((iR-1)*N+n,(iT-1)*N+n) = H_ls_i(n,iT,iR); 
                end
%                Heq((iR-1)*N+(1:N),(iT-1)*N+(1:N)) = diag(H_ls_i(:,iT,iR));    % SLOW
                
            end
            
        end     
        
        Hhat = H_ls_i;
    
    case 'LMMSE'  
        y = zeros(N,nR);
        Yp = zeros(p.Kon/p.Delta_k,nT,nR);
        H_lmmse_i = zeros(N,nT,nR);
        Heq = zeros(N*nT,N*nR);
        
        % noise
        VarS = 1; % Signal power
        VarN = VarS/(10^(snr/10));  
        R_NN = VarN*eye(p.Kon/p.Delta_k); 
        
        for iR = 1:nR
            % Receive signal in time
            y(:,iR) = ifft_u(Y(:,iR));
            %y(:,iR) = awgn(ifft_u(Y(:,iR)),snr);            
            for iT = 1:nT
                Yp(:,iT,iR) = Fp(:,:,iT)*y(:,iR);
                
                dXp = diag(Xp_iT(:,iT));
                
                R_YY = dXp*R_HH*dXp' + R_NN;
                
                H_lmmse = R_HH*dXp'/R_YY*Yp(:,iT,iR);   
                
                % FFT based interpolation through zero padding
                h_lmmse_temp = p.Delta_k*p.M*Fp(:,:,iT)'*H_lmmse/sqrt(N);                
                h_lmmse_temp(numPaths+1:end) = 0; 
                H_lmmse_i(:,iT,iR) = fft(h_lmmse_temp); 
                
                % vector to matrix
                for n = 1:N
                    Heq((iR-1)*N+n,(iT-1)*N+n) = H_lmmse_i(n,iT,iR); 
                end
%                Heq((iR-1)*N+(1:N),(iT-1)*N+(1:N)) = diag(H_lmmse_i(:,iT,iR)); % SLOW                
                
            end
        end
        
        Hhat = H_lmmse_i;         
   
end

