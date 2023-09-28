% Copyright (c) 2014 TU Dresden
% All rights reserved.
% See accompanying license.txt for details.
%
function h = get_mimo_channel(P,Fs,numPaths,path_delays,nT,nR,chan_type,fd)

% P: Power Delay Profile
% Fs: Sampling Frequency
% numPaths: Number of multipaths (number of channel taps)
% nT: Number of transmit antennas
% nR: Number of receive antennas
% fd: Doppler shift

if nargin < 7
    chan_type = 'block fading';
end

switch chan_type
    case 'block fading'
        g_TR = zeros(numPaths,nT,nR);
        h_i = zeros(numPaths,nT,nR);
        h = zeros(numPaths*nT,nR);
        for iR = 1:nR
            for iT = 1:nT
                g_TR(:,iT,iR) = complex(randn(numPaths, 1), randn(numPaths, 1))/sqrt(2);
                h_i(:,iT,iR) = sqrt(P(:,iT,iR)).*g_TR(:,iT,iR);
            end
            h(:,iR) = reshape(h_i(:,:,iR),[numPaths*nT,1]);
        end
    case 'time variant'
     %   Fs = 7.68e6; % Sampling frequency % 1.92 | 3.84 | 7.68 | 15.36 | 23.04 | 30.72
        %fd = 300;
        k =  0;  % K=0 equals rayleigh

        P_db = 10*log10(P(:,1,1));
        tau =  (path_delays-1)./Fs; % Path delays  (0:(length(P_db)-1))./Fs;
        
        clear h      
        for iR = 1:nR
            for iT = 1:nT                
                 chan_obj = ricianchan(1/Fs,fd,k,tau,P_db.');                 
                 chan_obj.ResetBeforeFiltering = 0;
              %   chan_obj.MaxDopplerShift = fd;
                 chan_obj.StoreHistory = 1;
                 chan_obj.StorePathGains = 1;
                 h(iT,iR) = chan_obj;
               %  warning('off','last')
            end           
        end
         
end

end
