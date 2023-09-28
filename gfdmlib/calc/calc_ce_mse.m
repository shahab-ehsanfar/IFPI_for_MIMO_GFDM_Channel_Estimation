function [ aEi, aEn, mEi, mEn, trD] = calc_ce_mse(p, Dd, Dp, H, Hhat, U, HFA, numPaths, snr, method, R_HH, R_XdXd, R_HHh, snrh, Fp, Delta_k, Y)
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

N = p.M*p.K;
if nargin < 10
   method = 'OP'; 
end

switch method
        %% One Pilot
    case 'OP'        
        % Initialization
        xp =  do_modulate(p, [Dp(:,1:p.M-1) Dp(:,p.M-1)]);
       
        Xp = fft(xp);
        
        Xp1 = 1./(Xp);
        Xpm = repmat(conj(Xp1), [1, p.K*p.M]);
        
        Z = Delta_k^2*(U'*U); 
        D = Xpm.*Z .* Xpm';   % Equal to: D = diag(Xp1)'*(Z)*diag(Xp1);
                
        % E[ Xd^H H^H D H Xd ]              
        sigma2d = diag(abs(Dd(:)) > 0) + 0;
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
        Dp1 = [Dp(:,1:p.M-1) Dp(:,p.M-1)];
        Dp2 = [Dp(:,2) Dp(:,2:p.M)];
        
        dd=Dd(:); dp1=Dp1(:); dp2=Dp2(:);
        
        Xd = fft(do_modulate(p,Dd));
        
        Xp1 =  fft(do_modulate(p, Dp1));
        Xp2 =  fft(do_modulate(p, Dp2));        
        
        Xp1inv = 1./(Xp1);
        Xp2inv = 1./(Xp2);
        
        Z = Delta_k^2*(U'*U); 
        
        D1 = diag(Xp1inv)'*(Z)*diag(Xp1inv);
        D2 = diag(Xp2inv)'*(Z)*diag(Xp2inv);
        D12 = diag(Xp1inv)'*(Z)*diag(Xp2inv);
              
            
        sigma2d = diag(abs(dd) > 0) + 0;
        % sigma2dp1 = diag(abs(dp1) > 0) + 0;
        % sigma2dp2 = diag(abs(dp2) > 0) + 0;
        
        Edp1d = (dp1*dp1') + sigma2d;
        Edp2d = (dp2*dp2') + sigma2d;
        
        %Edp1d = sigma2dp1 + sigma2d;
        %Edp2d = sigma2dp2 + sigma2d;
        
        SigmaPsi1Psi1 = HFA*Edp1d*HFA';
        SigmaPsi2Psi2 = HFA*Edp2d*HFA';
        aEi1 = trace((D1)*(SigmaPsi2Psi2));
        aEi2 = trace((D2)*(SigmaPsi1Psi1));
        
        SigmaPsi1Psi2 = HFA*( cov(dp2*dp1')+sigma2d )*HFA';
        
        D12s = 0.5*(D12+D12');
        Ei1i2 = trace(D12s*SigmaPsi1Psi2);
        
        aEi = real(0.25*(aEi1+aEi2 + 2*real(Ei1i2)));

        
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
        
    case 'fading'      
        % Initialization
        FA = HFA; % correct notation
        
        Xp = fft(do_modulate(p,Dp));
        
        Xp1 = 1./(Xp);
        Xpm = repmat(conj(Xp1), [1, p.K*p.M]);
        
        Z = Delta_k^2*(U'*U); 
        
        D = Xpm.*Z .* Xpm';   % Equal to: D = diag(Xp1)'*(Z)*diag(Xp1);
                   
%        sigma2d = diag(abs(Dd(:)) > 0) + 0;        
%         for n= 0:p.M*p.K-1
%            %hpath(n+1,:) = chan.PathGains.*sinc(n-tau); 
%            hpath(n+1,:) = sinc(n-tau);
%         end
%         P = diag(Gain.*sum(hpath,2));
%         
%         F = get_DFT_matrix(p);
%         FP = F*P*sqrt(p.M*p.K);
%         R_HH = FP*FP';        
        
%        R_XdXd = FA*sigma2d*FA';
                
        covPsi = R_XdXd .* R_HH;
        
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
        
        case 'LMMSE'
        % Initialization
        FA = HFA; % correct notation
        
        Xp = fft(do_modulate(p,Dp));
        
%         % E[ Xd^H H^H D H Xd ]
         sigma2d = diag(abs(Dd(:)) > 0) + 0;        
%         for n= 0:p.M*p.K-1
%            hpath(n+1,:) = sinc(n-tau);
%         end
%         P = diag(Gain.*sum(hpath,2));
%         
%         F = get_DFT_matrix(p);
%         FP = F*P*sqrt(p.M*p.K);
%         R_HH = FP*FP';              
        
%        R_XdXd = FA*sigma2d*FA';
                
        covPsi = R_XdXd .* R_HH;
        
        aEi = 0;
        
        % noise
        VarS = 1; % Signal power
        VarN = VarS/(10^(snr/10));
        
        R_NhNh = covPsi+VarN*eye(p.M*p.K);
        
        dXp = diag(Xp);
        %aEn= inv((dXp'/R_NhNh)*dXp+inv(R_HH));        
        R_YYinv = dXp*R_HH*dXp' + R_NhNh; %R_HH.*R_XdXd + VarN*eye(p.M*p.K); 
        aEn = trace(R_HH - R_HH*dXp'/R_YYinv*dXp*R_HH);
        
        mEn = sum(abs(Hhat - H).^2)./N;
        
        mEi = 0; 
        
        case 'LMMSEapprox'
        % Initialization
       % U = Z; % correct notation
 
        Xp = fft(do_modulate(p,Dp));        
       
        aEi = 0;
        
        % noise
        VarS = 1; % Signal power
        VarN = VarS/(10^(snr/10));
        VarNh = VarS/(10^(snrh/10));
        
        R_NhNh = R_XdXd .* R_HH +VarN*eye(p.M*p.K);
        R_NhNhp = R_XdXd .* R_HHh +VarNh*eye(p.M*p.K);
        
        dXp = diag(Xp);
           
        R_YY = dXp*R_HH*dXp' + R_NhNh; 
        R_YYh = dXp*R_HHh*dXp' + R_NhNhp; 
         
        w = R_HHh*dXp'/R_YYh;
               
        aEn = real(trace((U*w*R_YY*w'*U')+ R_HH - U*w*dXp*R_HH - R_HH*dXp'*w'*U'));
            
        mEn = sum(abs(Hhat - H).^2)./N;
        
        mEi = 0; 
        
        case 'LMMSEapproxNoL'
        % Initialization
         
        Xp = fft(do_modulate(p,Dp));        
       
        aEi = 0;
        
        % noise
        VarS = 1; % Signal power
        VarN = VarS/(10^(snr/10));
        VarNh = VarS/(10^(snrh/10));
        
        R_NhNh = R_XdXd .* R_HH +VarN*eye(p.M*p.K);
        R_NhNhp = R_XdXd .* R_HHh +VarNh*eye(p.M*p.K);
        
        dXp = diag(Xp);
           
        R_YY = dXp*R_HH*dXp' + R_NhNh; 
        R_YYh = dXp*R_HHh*dXp' + R_NhNhp; 
         
        w = R_HHh*dXp'/R_YYh;
               
        aEn = real(trace((w*R_YY*w')+ R_HH - w*dXp*R_HH - R_HH*dXp'*w'));
            
        mEn = sum(abs(Hhat - H).^2)./N;
        
        mEi = 0; 
        
        case 'OP_Delta_k'        
        % Initialization
%         Nf = (p.M*p.K)/Delta_k;
%                 
        H0 = H;
        H = Fp*ifft(H0)*sqrt(N);
                     
        %xp =  do_modulate(p, [Dp(:,1:p.M-1) Dp(:,p.M-1)]);
        xp = do_modulate(p, Dp);
        
        Xp = Fp*xp;
             
        Xp1 = 1./(Xp);
        Xpm = repmat(conj(Xp1), [1, p.K*p.M/Delta_k]);
        
%        Z = Delta_k^2*(U'*U); 
        
%        D = (Xpm.*Z .* Xpm');   % Equal to: D = diag(Xp1)'*(Z)*diag(Xp1);  
        
        D = U*diag(Xp1.*conj(Xp1))*U';

        % E[ Xd^H H^H D H Xd ]   
        
%         sigma2d = diag(abs(Dd(:)) > 0) + 0;
%         covPsi = HFA*sigma2d*HFA';
        covPsi = R_XdXd .* R_HH;
        
        aEi = real(trace((D)*covPsi))/(N/Delta_k);        
        
        % noise
        VarS = 1; % Signal power
        VarN = VarS/(10^(snr/10)); %VarS/(10^(snr/10)/Delta_k);
        aEn= real(VarN*trace(D))/(N/Delta_k);  
%        S = Fp*(VarN*eye(N))*Fp'/Delta_k;
%        aEn = real(trace(D*S))/(N/Delta_k);
        
        % Measured MSE
        xp0 =  do_modulate(p, Dp);
        xd0 = do_modulate(p, Dd);
         
        InterF = H.*((Fp*xd0)./(Fp*xp0));
        interf = Delta_k*Fp'*(InterF)/sqrt(N);   % Nf: makes Fp*Fp' unitary
        interf(numPaths+1:end)=0;
        InterFh = Fp*(interf)*sqrt(N);
        
        noise = Hhat - InterFh - H;
        mEn = sum(abs(noise).^2)./(N/Delta_k);
        
        mEi = sum(abs(InterFh).^2)./(N/Delta_k);
        
%         F = get_DFT_matrix(p);
%         Ft = Fp;
%         for i = numPaths+1:p.M*p.K
%             Ft(:,i) = zeros(p.M*p.K/Delta_k,1);
%         end
%         sum(abs(Delta_k*F*Ft'*Hhat - H0).^2)./length(H0)
%         sum(abs(Hhat - H).^2)./length(H)
%          []
        
    case 'LMMSE_Delta_k'          
        H0 = H;
        H = Fp*ifft(H0)*sqrt(N);
        
        xp = do_modulate(p, Dp);
        
        Xp = Fp*xp;          
       
        %R_HH = R_HH/(N/Delta_k);
        covPsi = R_XdXd .* R_HH;
        
        aEi = 0;  
        
        % noise
%         VarS = 1; % Signal power
%         VarN = VarS/(10^(snr/10));
        VarN = 10^(-snr/10);
        S = VarN*eye(N/Delta_k); %Fp*(VarN*eye(N))*Fp'/(N/Delta_k); % N/Delta_k normalizes ifft
        
        R_NhNh = covPsi + S;
        
        dXp = diag(Xp);   
        %aEn= inv((dXp'/R_NhNh)*dXp+inv(R_HH));        
        R_YYinv = dXp*R_HH*dXp' + R_NhNh; %R_HH.*R_XdXd + VarN*eye(p.M*p.K); 
        aEn = trace(R_HH - R_HH*dXp'/R_YYinv*dXp*R_HH);
        aEn = real(aEn)/(N/Delta_k);
        
        mEn = sum(abs(Hhat - H).^2)./(N/Delta_k);
        
        mEi = 0; 
        
    case 'HM_Delta_k'
        %% Harmonic Mean
        % Initialization
        FA = HFA; %correct notation
        H0 = H;
        H = Fp*ifft(H0)*sqrt(N);
        
        Dp1 = [Dp(:,1:p.M-1) Dp(:,p.M-1)];
        Dp2 = [Dp(:,2) Dp(:,2:p.M)];
        
        dd=Dd(:); dp1=Dp1(:); dp2=Dp2(:);
        
        Xd = Fp*(do_modulate(p,Dd));
        
        Xp1 =  Fp*(do_modulate(p, Dp1));
        Xp2 =  Fp*(do_modulate(p, Dp2));
        
        Xp1inv = 1./(Xp1);
        Xp2inv = 1./(Xp2);
        
        Z = Delta_k^2*(U'*U); 
        
        D1 = diag(Xp1inv)'*(Z)*diag(Xp1inv);
        D2 = diag(Xp2inv)'*(Z)*diag(Xp2inv);
        D12 = diag(Xp1inv)'*(Z)*diag(Xp2inv);
        
       
        sigma2d = diag(abs(dd) > 0) + 0;
        % sigma2dp1 = diag(abs(dp1) > 0) + 0;
        % sigma2dp2 = diag(abs(dp2) > 0) + 0;
        
        Edp1d = (dp1*dp1') + sigma2d;
        Edp2d = (dp2*dp2') + sigma2d;
        
        Edp1d = sigma2d;
        Edp2d = sigma2d;
        
        %Edp1d = sigma2dp1 + sigma2d;
        %Edp2d = sigma2dp2 + sigma2d;
        
        SigmaX1X1 = FA*Edp1d*FA';
        SigmaX2X2 = FA*Edp2d*FA';
        
        SigmaPsi1Psi1 = SigmaX1X1 .* R_HH;
        SigmaPsi2Psi2 = SigmaX2X2 .* R_HH;
        aEi1 = trace((D1)*(SigmaPsi2Psi2)) + trace(diag(Xp2)'*D1*diag(Xp2)*R_HH);
        aEi2 = trace((D2)*(SigmaPsi1Psi1)) + trace(diag(Xp1)'*D2*diag(Xp1)*R_HH);
        
        %SigmaX1X2 = FA*( cov(dp2*dp1')+sigma2d )*FA';
        SigmaX1X2 = FA*(sigma2d)*FA';
        SigmaPsi1Psi2 = SigmaX1X2 .* R_HH;
                
        %D12s = 0.5*(D12+D12');
        Ei1i2 = trace(D12*SigmaPsi1Psi2) + trace(diag(Xp2)'*D12*diag(Xp1)*R_HH);
        
        aEi = real(0.25*(aEi1+aEi2 + 2*real(Ei1i2)))/(N/Delta_k);
        
        
        % noise
        VarS = 1; % Signal power
        VarN = VarS/(10^(snr/10)); % equivalent to S in 'OP_Delta_k'
                
        aEn1= real(VarN*trace(D1));
        aEn2= real(VarN*trace(D2));
        aEn12 = (VarN*trace(D12));
        
        aEn = 0.25*(aEn1+aEn2 + 2*real(aEn12))/(N/Delta_k);
        
        % Measured MSE
        InterF = 0.5*H.*(Xd+Xp2)./Xp1 + 0.5*H.*(Xd+Xp1)./Xp2;
        interf = Delta_k*Fp'*(InterF)/(sqrt(N));
        interf(numPaths+1:end)=0;
        InterFh = Fp*(interf)*sqrt(N);
        
        noise = Hhat - InterFh - H;
        mEn = sum(abs(noise).^2)./(N/Delta_k);
        mEi = sum(abs(InterFh).^2)./(N/Delta_k);

     %   Xp = Xp1 + Xp2;
     %   trD = abs([trace((diag(1./Xp)'*Z*diag(1./Xp))) trace(D1) trace(0.25*(D1+D2+2*D12))]).';
end

end




% function [ aEi, aEn, mEi, mEn] = calc_ce_mse(p, Dd, Dp, H, Hhat, Z, HFA, numPaths, snr, method)
% %
% % Analytical MSE computation
% % Dd, Dp are data and pilot matrices respectively
% % 
% % FA: DFT Matrix times transmitter matrix A
% % Hhat: Estimated channel frequency response
% % H: Channel frequency response
% % Z: (F_t F^H)^H (F_t F^H) = U^H U
% %
% % aEi: Analytic MSE due to interference  
% % aEn: Analytic MSE due to noise
% % mEi: Measured MSE due to interference
% % mEn: Measured MSE due to noise
% % 
% % Author: Shahab Ehsanfar
% 
% if nargin < 10
%    method = 'OP'; 
% end
% 
% switch method
%         %% One Pilot
%     case 'OP'         
%         % Initialization
%         Dd(:,1) = zeros(p.K,1);
%         
%         N = p.M*p.K;
%         
%         Xd = fft(do_modulate(p,Dd));
%         Xp = fft(do_modulate(p,Dp));
%         
%         Xp1 = 1./(Xp);
%         Xpm = repmat(conj(Xp1), [1, p.K*p.M]);
%         D = Xpm.*Z .* Xpm';   % Equal to: D = diag(Xp1)'*(Z)*diag(Xp1);
%                 
%         % E[ Xd^H H^H D H Xd ]        
%         %sigma2d = diag(var(Dd(:)*Dd(:)'));        
%         sigma2d = diag(abs(Dd(:)) > 0) + 0;
%         tic;
%         covPsi = HFA*sigma2d*HFA';        
%         toc
%         aEi = real(trace((D)*covPsi));        
%         
%         % noise
%         VarS = 1; % Signal power
%         VarN = VarS/(10^(snr/10));
%         aEn= real(VarN*trace(D));        
%         
%         % Measured MSE
%         xp0 =  do_modulate(p, Dp);
%         xd0 = do_modulate(p, Dd);
%         InterF = H.*(fft(xd0)./fft(xp0));
%         interf = ifft(InterF);
%         interf(numPaths+1:end)=0;
%         InterFh = fft(interf);
%         
%         noise = Hhat - InterFh - H;
%         mEn = sum(abs(noise).^2)./N;
%         
%         mEi = sum(abs(InterFh).^2)./N;
% 
%         
%     case 'HM' 
%         %% Harmonic Mean
%         % Initialization
%         Dd(:,1) = zeros(p.K,1);
%         Dd(:,end) = zeros(p.K,1);
%         Dp1 = [Dp(:,1:p.M-1) Dp(:,p.M-1)];
%         Dp2 = [Dp(:,2) Dp(:,2:p.M)];
%         
%         dd=Dd(:); dp1=Dp1(:); dp2=Dp2(:);
%         
%         N = p.M*p.K;
%         
%         Xd = fft(do_modulate(p,Dd));
%         
%         Xp1 =  fft(do_modulate(p, Dp1));
%         Xp2 =  fft(do_modulate(p, Dp2));        
%         
%         Xp1inv = 1./(Xp1);
%         Xp2inv = 1./(Xp2);
%         
%         D1 = diag(Xp1inv)'*(Z)*diag(Xp1inv);
%         D2 = diag(Xp2inv)'*(Z)*diag(Xp2inv);
%         D12 = diag(Xp1inv)'*(Z)*diag(Xp2inv);
%               
%             
%         sigma2d = diag(abs(dd) > 0) + 0;
%         
%         Edp1d = cov(dp1*dp1') + sigma2d;
%         Edp2d = cov(dp2*dp2') + sigma2d;
%         
%         SigmaPsi1Psi1 = HFA*Edp1d*HFA';
%         SigmaPsi2Psi2 = HFA*Edp2d*HFA';
%         aEi1 = trace((D1)*(SigmaPsi2Psi2));
%         aEi2 = trace((D2)*(SigmaPsi1Psi1));
%         
%         SigmaPsi1Psi2 = HFA*( cov(dp2*dp1')+sigma2d )*HFA';
%         
%         D12s = 0.5*(D12+D12');
%         Ei1i2 = trace(D12s*SigmaPsi1Psi2);
%         
%         aEi = real(0.25*(aEi1+aEi2 + 2*real(Ei1i2)));
% 
%         
%         % noise
%         VarS = 1; % Signal power
%         VarN = VarS/(10^(snr/10));
%         aEn1= real(VarN*trace(D1));
%         aEn2= real(VarN*trace(D2));
%         aEn12 = (VarN*trace(D12s));
%         
%         aEn = 0.25*(aEn1+aEn2 + 2*real(aEn12));
%         
%         % Measured MSE
%         InterF = 0.5*H.*(Xd+Xp2)./Xp1 + 0.5*H.*(Xd+Xp1)./Xp2;
%         interf = ifft(InterF);
%         interf(numPaths+1:end)=0;
%         InterFh = fft(interf);
%         
%         noise = Hhat - InterFh - H;
%         mEn = sum(abs(noise).^2)./N;
%         mEi = sum(abs(InterFh).^2)./N;
%         
% end
% 
% end
% 

