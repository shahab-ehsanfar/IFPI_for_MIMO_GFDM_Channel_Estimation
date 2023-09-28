function [ Hhat ] = do_Chan_Estimation( p, Y, Dp, numPaths, method, snr, R_HH, R_XdXd, H, Dd, Fp, Delta_k)
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


% Separate the first and second pilot subsymbols

if p.M > 1
    xp1 =  do_modulate(p, [Dp(:,1:p.M-1) Dp(:,p.M-1)]);    
else 
    xp1 =  do_modulate(p, Dp);    
end

switch method
    case 1
        % Estimate the channel through LS using only one pilot subsymbol
     %   xp1 = do_modulate(p, Dp);
        Xp1 = fft(xp1);
        H_ls = (Y./Xp1);
        
        % Apply the a-priori-Knowledge of number of paths to the estimation 
        h_ls = ifft(H_ls);
        h_ls(numPaths+1:end) = 0;
        Hhat = fft(h_ls); 
        
    case 2        
        % Estimate the channel through Harmonic mean of the pilot subsymbols
        Xp1 = fft(xp1);
        xp2 =  do_modulate(p, [Dp(:,2) Dp(:,2:p.M)]);        
        Xp2 = fft(xp2);
        
        H_ls = 0.5*Y.*(1./Xp1 + 1./Xp2);
        
        % Apply the a-priori-Knowledge of number of paths to the estimation (1)
        h_ls = ifft(H_ls);
        h_ls(numPaths+1:end) = 0;
        Hhat = fft(h_ls);
    case 3
        % LMMSE estimation for one pilot subsymbol
        Xp1 = fft(xp1);
        R_PsiPsi = R_HH.*R_XdXd;
       
        % noise
        VarS = 1; % Signal power
        VarN = VarS/(10^(snr/10));                      
        
        dXp = diag(Xp1);
        R_NN = R_PsiPsi + VarN*eye(p.M*p.K);       
                
        invFactor = dXp*R_HH*dXp' + R_NN;
        H_lmmse = R_HH*dXp'/invFactor*Y;        
        
        % Apply the a-priori-Knowledge of number of paths to the estimation
        h_lmmse = ifft(H_lmmse);
        h_lmmse(numPaths+1:end) = 0;
        Hhat = fft(h_lmmse); % the final estimation
    case 4
        % TWO pilots LMMSE
        % Estimation for First-Last One Pilot arrangement        
        Xp = fft(do_modulate(p, Dp));
        
        R_PsiPsi = R_HH.*R_XdXd;
       
        % noise
        VarS = 1; % Signal power
        VarN = VarS/(10^(snr/10));                      
        
        dXp = diag(Xp);
        R_NN = R_PsiPsi + VarN*eye(p.M*p.K);       
                
        invFactor = dXp*R_HH*dXp' + R_NN;
        H_lmmse = R_HH*dXp'/invFactor*Y;        
        
        % Apply the a-priori-Knowledge of number of paths to the estimation
        h_lmmse = ifft(H_lmmse);
        h_lmmse(numPaths+1:end) = 0;
        Hhat = fft(h_lmmse); % the final estimation
    case 5
        % LMMSE estimation whithout a-priory-Knowledge
        Xp1 = fft(xp1);
        R_PsiPsi = R_HH.*R_XdXd;
       
        % noise
        VarS = 1; % Signal power
        VarN = VarS/(10^(snr/10));                      
        
        dXp = diag(Xp1);
        R_NN = R_PsiPsi + VarN*eye(p.M*p.K);       
                
        invFactor = dXp*R_HH*dXp' + R_NN;
        Hhat = R_HH*dXp'/invFactor*Y;        
        
    case 6
        % CP: HM estimation by removing CP part
        xp2 =  do_modulate(p, [Dp(:,2) Dp(:,2:p.M)]);                      
        Xp1_cp1 = fft(xp1); Xp2_cp1 = fft(xp2);
               
        y = ifft(Y);
        Y1 = fft(y(p.K+1:p.M*p.K+p.K)); 
        
        H_ls1 = 0.5*Y1.*(1./Xp1_cp1 + 1./Xp2_cp1);
                
        h_ls1 = ifft(H_ls1);
        h_ls1(numPaths+1:end) = 0;
        Hhat = fft(h_ls1(1:end),p.K*p.M+p.K);
        
    case 7
        % Double Harmonic Mean
        % CP: HM estimation by considering also CP part 
        xp2 =  do_modulate(p, [Dp(:,2) Dp(:,2:p.M)]);  
        xp1_cp = [xp1(p.M*p.K-p.K+1:p.M*p.K).',xp1.'].';
        xp2_cp = [xp2(p.M*p.K-p.K+1:p.M*p.K).',xp2.'].';
        
        xp1_cp1 = xp1; xp2_cp1 = xp2;
        Xp1_cp1 = fft(xp1_cp1); Xp2_cp1 = fft(xp2_cp1);
        
        xp1_cp2 = xp1_cp(1:p.K*p.M); xp2_cp2 = xp2_cp(1:p.K*p.M);
        Xp1_cp2 = fft(xp1_cp2); Xp2_cp2 = fft(xp2_cp2);
        
        y = ifft(Y);
        Y1 = fft(y(p.K+1:p.K*p.M+p.K)); Y2 = fft(y(1:p.M*p.K));
               
        H_ls_mean = 0.25*(Y1.*(1./Xp1_cp1 + 1./Xp2_cp1)+Y2.*(1./Xp1_cp2 + 1./Xp2_cp2));
        
        h_ls_m = ifft(H_ls_mean);
        h_ls_m(numPaths+1:end) = 0;
        Hhat = fft(h_ls_m(1:end),p.K*p.M+p.K);   
       
    case 8
        % CP: consider only the edges
        xp2 =  do_modulate(p, [Dp(:,2) Dp(:,2:p.M)]);  
        xp1_cp = [xp1(p.M*p.K-p.K+1:p.M*p.K).',xp1.'].';
        xp2_cp = [xp2(p.M*p.K-p.K+1:p.M*p.K).',xp2.'].';
               
        xp_s = [xp1_cp(p.K*p.M+1:p.K*p.M+p.K);xp2_cp(1:p.K)];
        Xp_s1 = fft(xp_s);
        
        xp_s2 = [xp2_cp(p.K*p.M+1:p.K*p.M+p.K);xp1_cp(1:p.K)];
        Xp_s2 = fft(xp_s2);           

        y = ifft(Y);
        y_0 = [y(p.K*p.M+1:p.K*p.M+p.K);y(1:p.K)];
        Y_est = fft(y_0);       
                       
        H_ls_new = 0.5*Y_est.*(1./Xp_s1 + 1./Xp_s2);
        
        h_ls_new = ifft(H_ls_new);
        h_ls_new(numPaths+1:end) = 0;
        Hhat = fft(h_ls_new(1:end),p.K*p.M+p.K);
     case 9
        % Cyclic Prefix LMMSE
        xp = (do_modulate(p, Dp));
        xp_cp = [xp(p.M*p.K-p.K+1:p.M*p.K);xp];
        Xp = fft(xp_cp);
        
        R_PsiPsi = R_HH.*R_XdXd;
       
        % noise
        VarS = 1; % Signal power
        VarN = VarS/(10^(snr/10));                      
        
        dXp = diag(Xp);
        R_NN = R_PsiPsi + VarN*eye(p.M*p.K+p.K);       
                
        invFactor = dXp*R_HH*dXp' + R_NN;
        H_lmmse = R_HH*dXp'/invFactor*Y;        
        
%         % Apply the a-priori-Knowledge of number of paths to the estimation
%         h_lmmse = ifft(H_lmmse);
%         h_lmmse(numPaths+1:end) = 0;
%         Hhat = fft(h_lmmse); % the final estimation
        Hhat = H_lmmse;
        
    case 10
        % Estimate the channel through LS using only one pilot subsymbol
%         H0 = H;
%         H = Fp*ifft(H0);
        xp1 = do_modulate(p, Dp);    
        N = p.M*p.K;
        Xp1 = Fp*xp1;
        H_ls = (Y./Xp1);

        
        % Apply the a-priori-Knowledge of number of paths to the estimation
        h_ls = Delta_k*Fp'*(H_ls)/sqrt(N); 
        h_ls(numPaths+1:end) = 0;
        %Hhat = Fp*(h_ls)*sqrt(N); 
        Hhat = h_ls;   
        
    case 11
        N = p.M*p.K;
        xp = do_modulate(p, Dp);        
        Xp = Fp*xp;      
        
        R_PsiPsi = R_HH.*R_XdXd;
       
        % noise
        VarS = 1; % Signal power
        VarN = VarS/(10^(snr/10));  
        S = VarN*eye(N/Delta_k); %Fp*(VarN*eye(N))*Fp'/(N/Delta_k);
        
        dXp = diag(Xp);
        R_NN = R_PsiPsi + S;       
                
        invFactor = dXp*R_HH*dXp' + R_NN;
        H_lmmse = R_HH*dXp'/invFactor*Y;   
        
        if 0 % TEST
        H_lmmse2 = R_HH/(R_HH + R_NN/(dXp'*dXp))*(dXp\Y);
        figure; plot(abs(H_lmmse2))
        H_lmmse3 = R_HH*inv(R_HH+R_NN*inv(dXp'*dXp))*(inv(dXp)*Y);
        figure; plot(abs(H_lmmse3))
        end
       
        % Apply the a-priori-Knowledge of number of paths to the estimation
        h_lmmse = Delta_k*Fp'*(H_lmmse)/sqrt(N);
        h_lmmse(numPaths+1:end) = 0;
%         Hhat = Fp*(h_lmmse); % the final estimation
        %Hhat = H_lmmse;
        Hhat = h_lmmse;
        
        %Hp = Fp*ifft(H);
        
    case 12
        N = p.M*p.K;
        Xp1 = Fp*(xp1);
        xp2 =  do_modulate(p, [Dp(:,2) Dp(:,2:p.M)]);        
        Xp2 = Fp*(xp2);
        
        H_hm = 0.5*Y.*(1./Xp1 + 1./Xp2);
        
        % Apply the a-priori-Knowledge of number of paths to the estimation (1)
        h_hm = Delta_k*Fp'*(H_hm)/(sqrt(N));
        h_hm(numPaths+1:end) = 0;
        %Hhat = Fp*(h_ls)*sqrt(N);
        Hhat = h_hm;
                           
end

end

