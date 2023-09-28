function corr = ref_sync_autocorrelate(distance, rx_signal, Ncp, Ncs, mode, type, window)
% Perform the autocorrelation of the signal with itself
% but only calculate autocorrelation between config.preamble.K
% adjacent samples

if nargin <6 
   type = 'default';
end

if strcmp(type,'default')    
    K = distance;
%     shift_1 = rx_signal(1:end-K);
%     shift_2 = rx_signal((K+1):end);  
    
    
    shift_1 = [zeros(K,1); rx_signal];
    shift_2 = [rx_signal; zeros(K,1)];
    shift_2_mag = [rx_signal; 1*rx_signal(1:K)];
 
    
    if (nargin < 5 || strcmp(mode,'AC'))
        mult = conj(shift_1) .* (shift_2);
    elseif strcmp(mode,'CA') % Conjugate Autocorrelation
        mult = shift_1 .* (shift_2);
    elseif strcmp(mode,'RA') % Reverse Autocorrelation
        mult = shift_1 .* (flipud(shift_2));
    end
    
    autocorrelation = (conv(mult, ones(K,1)));
    magnitude = abs(conv(abs(shift_2_mag).^2, ones(K,1))).^1; 
    
    corr = autocorrelation.^1 ./ (magnitude+1).^1;
    
    if Ncp == 0
        % new line for CP-less
%         corr = circshift(corr,-K);
        corr = corr(K+K:end-K);
    end
    
    if (Ncp+Ncs > 0)
        CpCs = Ncp + Ncs;
        corr = conv(corr, ones(CpCs,1)) / CpCs;
        
        corr = corr((2*K-CpCs):end);
    end
elseif strcmp(type,'manual')
    K = distance;
    W = window; % Window size
    N_rx_signal = length(rx_signal);
             
    autocorr = zeros(N_rx_signal,1);
    magnitude = zeros(N_rx_signal,1);
    for n = 1:N_rx_signal-K-W        
        autocorr(n) = sum( conj(rx_signal(n:n+W-1)).*rx_signal(n+K:n+K+W-1) );
        magnitude(n) = sum( abs(rx_signal(n+K:n+K+W-1)).^2 );
    end    
    
    corr = (abs(autocorr).^2)./(abs(magnitude+1).^2);
elseif strcmp(type,'different_window')
    K = distance;
       
    shift_1 = [zeros(K,1); rx_signal];
    shift_2 = [rx_signal; zeros(K,1)];
    shift_2_mag = [rx_signal; 0.4*rx_signal(1:K)];
        
    if (nargin < 5 || strcmp(mode,'AC'))
        mult = conj(shift_1) .* (shift_2);
    elseif strcmp(mode,'CA') % Conjugate Autocorrelation
        mult = shift_1 .* (shift_2);
    elseif strcmp(mode,'RA') % Reverse Autocorrelation
        mult = shift_1 .* (flipud(shift_2));
    end
    
    autocorrelation = (conv(mult, ones(window,1)));
    magnitude = abs(conv(abs(shift_2_mag).^2, ones(window,1))).^1; % abs(conv(abs(shift_2).^2, ones(K+49,1))).^2;
    
    corr = autocorrelation.^1 ./ (magnitude+1).^1;
    
    if Ncp == 0
        % new line for CP-less
        %corr = circshift(corr,-K);
        corr = corr(K+window:end-K);
    end
    
    if (Ncp+Ncs > 0)
        CpCs = Ncp + Ncs;
        corr = conv(corr, ones(CpCs,1)) / CpCs;
        
        corr = corr((2*K-CpCs):end);
    end
end

end

