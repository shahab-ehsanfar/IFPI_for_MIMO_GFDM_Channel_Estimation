function [ aEi, aEn, mEi, mEn] = calc_ce_mse(p, Dd, Dp, H, Hhat, Z, HFA, numPaths, snr, method)
%
% Analytical MSE computation
% Dd, Dp are data and pilot matrices respectively
% 
% FA: DFT Matrix times transmitter matrix A
% Hhat: Estimated channel frequency response
% H: Channel frequency response
% Z: (F_t F^H)^H (F_t F^H) = U^H U
%
% aEi: Analytic MSE due to interference  
% aEn: Analytic MSE due to noise
% mEi: Measured MSE due to interference
% mEn: Measured MSE due to noise
% 
% Author: Shahab Ehsanfar

if nargin < 10
   method = 'OP'; 
end

switch method
        %% One Pilot
    case 'OP'         
        % Initialization
        L=numPaths;
        Dd(:,1) = zeros(p.K,1);
        
        N = p.M*p.K;
        
        Xd = fft(do_modulate(p,Dd));
        Xp = fft(do_modulate(p,Dp));
        
        Xp1 = 1./(Xp);
        
        % E[ Xd^H H^H D H Xd ]
        D = diag(Xp1)'*(Z)*diag(Xp1);
        
        psi = diag(H)*Xd;
        
        %mEi = abs(psi'*(D)*psi/N);
        
        sigma2d = diag(var(Dd(:)*Dd(:)'));
        
        %covPsi = diag(H)*FA*sigma2d*FA'*diag(H)';
        covPsi = HFA*sigma2d*HFA';
        
        aEi = real(trace((D)*covPsi));
        
        
        % noise
        VarS = 1; % Signal power
        VarN = VarS/(10^(snr/10));
        aEn= real(VarN*trace(D));
        
        
        % Measured MSE
        xp0 =  do_modulate(p, Dp);
        xd0 = do_modulate(p, Dd);
        InterF = H.*(fft(xd0)./fft(xp0));
        interf = ifft(InterF);
        interf(numPaths+1:end)=0;
        InterFh = fft(interf);
        
        noise = Hhat - InterFh - H;
        mEn = sum(abs(noise).^2)./N;
        
        mEi = sum(abs(InterFh).^2)./N;

        
    case 'HM' 
        %% Harmonic Mean
        % Initialization
        L=numPaths;
        Dd(:,1) = zeros(p.K,1);
        Dd(:,end) = zeros(p.K,1);
        Dp1 = [Dp(:,1:p.M-1) Dp(:,p.M-1)];
        Dp2 = [Dp(:,2) Dp(:,2:p.M)];
        
        dd=Dd(:); dp1=Dp1(:); dp2=Dp2(:);
        
        N = p.M*p.K;
        
        Xd = fft(do_modulate(p,Dd));
        
        Xp1 =  fft(do_modulate(p, Dp1));
        Xp2 =  fft(do_modulate(p, Dp2));        
        
        Xp1inv = 1./(Xp1);
        Xp2inv = 1./(Xp2);
        
%         Xd1 = Xp2 + Xd;
%         Xd2 = Xp1 + Xd;
        
        %HFA = diag(H)*FA;
        
        % E[ Xd^H H^H D H Xd ]
        D1 = diag(Xp1inv)'*(Z)*diag(Xp1inv);
        D2 = diag(Xp2inv)'*(Z)*diag(Xp2inv);
        D12 = diag(Xp1inv)'*(Z)*diag(Xp2inv);
        
%         psi1 = diag(H)*Xd1;
%         psi2 = diag(H)*Xd2;
%         
%         mEi1b = (psi1'*(D1)*psi1);
%         mEi2b = (psi2'*(D2)*psi2);
%         
%         % INNER PRODUCT EXACT
%         I12inner = Xd1'*diag(H)'*diag(Xp1inv)'*Z*diag(Xp2inv)*diag(H)*Xd2;
%         Iinner = I12inner + I12inner';
%         
%         mEib = 0.25*abs(mEi1b + mEi2b + Iinner)/N;
%         % % % % % 
        
%         varddp = var(Dp2(:)*Dp1(:)');
%         vardd = var(Dd(:)*Dd(:)');
%         
%         %covPsi12 = diag(H)*FA*diag(varddp+vardd)*FA'*diag(H)'; 
%         covPsi12 = cov(psi1*psi2')/N^2;
%         Ei1i2 = trace((0.5*(D12+D12'))*covPsi12);

        SigmaPhiPhi = HFA*( cov(dp2*dp1')+diag(var(dd*dd')) )*HFA';
%        SigmaPsiPsi = HFA*diag(var(Dd(:)*Dd(:)'))*HFA';
        
        D12s = 0.5*(D12+D12');
        Ei1i2 = trace(D12s*SigmaPhiPhi); %+ trace(D12s*SigmaPsiPsi);
        
        
%         % THE FIRST TWO TERMS
%         sigma2d1 = diag(var(dd1*dd1'));
%         sigma2d2 = diag(var(dd2*dd2'));
%         
%         covPsi1 = diag(H)*FA*sigma2d1*FA'*diag(H)';
%         covPsi2 = diag(H)*FA*sigma2d2*FA'*diag(H)';
%         aEi10 = trace((D1)*covPsi1);
%         aEi20 = trace((D2)*covPsi2);
%         aEi = 0.25*(aEi10+aEi20 + 2*real(Ei1i2));
        
        Edp1d = diag( var(dp1*dp1') + var(dd*dd') );
        Edp2d = diag( var(dp2*dp2') + var(dd*dd') );
%         SigmaPhiPhi1 = HFA*diag(var(Dp1(:)*Dp1(:)'))*HFA';
%         SigmaPhiPhi2 = HFA*diag(var(Dp2(:)*Dp2(:)'))*HFA';
%         
% %         aEi1 = trace((D1)*SigmaPsiPsi)+ trace((D1)*SigmaPhiPhi2);
% %         aEi2 = trace((D2)*SigmaPsiPsi)+ trace((D2)*SigmaPhiPhi1);
%         aEi1 = trace((D1)*(SigmaPsiPsi+SigmaPhiPhi2));
%         aEi2 = trace((D2)*(SigmaPsiPsi+SigmaPhiPhi1));
        SigmaPhiPhi1 = HFA*Edp1d*HFA';
        SigmaPhiPhi2 = HFA*Edp2d*HFA';
        aEi1 = trace((D1)*(SigmaPhiPhi2));
        aEi2 = trace((D2)*(SigmaPhiPhi1));
%         
%         % aEi = abs(0.25*(aEi1+aEi2+Iinner/N));
%         %aEi = 0.25*(aEi1+aEi2 + Ei1i2 + Ei1i2');
        aEi = real(0.25*(aEi1+aEi2 + 2*real(Ei1i2)));
%         aEi = abs(real(aEi));
        
        % noise
        VarS = 1; % Signal power
        VarN = VarS/(10^(snr/10));
        aEn1= real(VarN*trace(D1));
        aEn2= real(VarN*trace(D2));
        aEn12 = (VarN*trace(D12s));
        
        aEn = 0.25*(aEn1+aEn2 + 2*real(aEn12));
        
        
        % Measured MSE
        InterF = 0.5*H.*(Xd+Xp2)./Xp1 + 0.5*H.*(Xd+Xp1)./Xp2;
        interf = ifft(InterF);
        interf(numPaths+1:end)=0;
        InterFh = fft(interf);
        
        noise = Hhat - InterFh - H;
        mEn = sum(abs(noise).^2)./N;
        mEi = sum(abs(InterFh).^2)./N;
        
end

end

