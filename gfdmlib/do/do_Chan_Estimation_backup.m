function [ Hhat ] = do_Chan_Estimation( p, Y, Dp, numPaths, method)
%
% method:
%       1: Estimate the channel using only one pilot subsymbol
%       2: Estimate the channel through Harmonic mean of the pilot subsymbols
%       3: Performs Channel Estimation using Maximum Likelihood approach
%
%
% X:    FFT of transmit signal
% Y:    FFT of received signal
% Dp:   K*M Matrix of Pilot signals with zero on data positions
% numPaths: Number of paths in FSC situation
%
% Author: Shahab Ehsanfar


% Separate the first and second pilot subsymbols
xp1 =  do_modulate(p, [Dp(:,1:p.M-1) Dp(:,p.M-1)]);
xp2 =  do_modulate(p, [Dp(:,2) Dp(:,2:p.M)]);

Xp1 = fft(xp1);
Xp2 = fft(xp2);


switch method
    case 1
        % Estimate the channel using only one pilot subsymbol
        H_ls = (Y./Xp1);
        
        % Apply the a-priori-Knowledge of number of paths to the estimation 
        h_ls = ifft(H_ls);
        h_ls(numPaths+1:end) = 0;
        Hhat = fft(h_ls); 
        
    case 2        
        % Estimate the channel through Harmonic mean of the pilot subsymbols
        H_ls = 0.5*Y.*(1./Xp1 + 1./Xp2);
        
        % Apply the a-priori-Knowledge of number of paths to the estimation (1)
        h_ls = ifft(H_ls);
        h_ls(numPaths+1:end) = 0;
        Hhat = fft(h_ls);    
    
    case 3 % Performs Channel Estimation using Maximum Likelihood approach
        %        
        % Estimate the channel through first pilot using LS (1)
        H_ls = Y./Xp1;
                
        % Apply the a-priori-Knowledge of number of paths to the estimation (1)
        h_ls = ifft(H_ls);
        h_ls(numPaths+1:end) = 0;
        Ha = fft(h_ls);
        
        % Estimate the channel through second pilot using LS (2)
        H_ls2 = Y./Xp2;
        
        % Apply the a-priori-Knowledge of number of paths to the estimation (2)
        h_ls2 = ifft(H_ls2);
        h_ls2(numPaths+1:end) = 0;
        Ha2 = fft(h_ls2);
        
        % Perform the Maximum Likelihood estimation
        Hhat= zeros(length(Y),1);
        Xhat= zeros(length(Y),1);
              
        for i= 1:length(Y)
            ML_1 = abs(Y(i) - Xp1(i).*Ha(i)).^2;
            ML_2 = abs(Y(i) - Xp2(i).*Ha2(i)).^2;
        
            if (ML_1 < ML_2) 
                Xhat(i) = Xp1(i);
            else
                Xhat(i) = Xp2(i);
            end
        end 
        
        Hhat = Y./Xhat;
        
        for i= 1:length(Y)
            ML_1 = abs(Y(i) - Xhat(i).*Ha(i)).^2;
            ML_2 = abs(Y(i) - Xhat(i).*Ha2(i)).^2;
        
            if (ML_1 < ML_2) 
                Hhat(i) = Ha(i);
            else
                Hhat(i) = Ha2(i);
            end
        end        
        
        
        
end

end

