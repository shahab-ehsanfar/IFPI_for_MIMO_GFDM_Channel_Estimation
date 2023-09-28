function [ Hhat ] = do_Chan_Estimation( p, Y, Dp, numPaths, method, hm)
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

if nargin == 6
    
    Dp1 = [Dp(:,1:p.M-1) Dp(:,p.M-1)];    
    Dp2 = [Dp(:,2) Dp(:,2:p.M)];
    
    xp1 =  do_modulate(p, Dp1);
    xp2 =  do_modulate(p, Dp2);
    
    Xp1 = fft(xp1,p.M*p.K);
    Xp2 = fft(xp2,p.M*p.K);
else
    xp1 =  do_modulate(p, [Dp(:,1:p.M-1) Dp(:,p.M-1)]);
    xp2 =  do_modulate(p, [Dp(:,2) Dp(:,2:p.M)]);
    
    Xp1 = fft(xp1);
    Xp2 = fft(xp2);
end


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
    case 3
        H_ls = (Y./Xp1);
        
        h_ls = ifft(H_ls);
        
        h_ls(numPaths+1:end) = 0;
        Hhat_0 = fft(h_ls);
        
        % Apply the a-priori-Knowledge of number of paths to the estimation
        
        HXdhat = Y - Hhat_0.* Xp1;
        Hhat1_2 = Hhat_0 - (HXdhat./Xp1);
        h_ls = ifft(Hhat1_2);
        
        h_ls(numPaths+1:end) = 0;
        Hhat = fft(h_ls);
        

% MSE_1 = sum(abs(H - Hhat1).^2)
% MSE_2 = sum(abs(H - Hhat1_2).^2)

        
end

end

