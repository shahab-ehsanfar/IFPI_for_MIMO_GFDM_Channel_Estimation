function crosscorr = ref_sync_crosscorrelation(preamble, rx_signal, type, N, noUW)

if nargin<3
   type = 'default_old';   
end

if strcmp(type, 'default_old')
    P = length(preamble);
    crosscorr = abs(conv(rx_signal, conj(preamble)) / sqrt((preamble'*preamble))).^2;
    crosscorr = crosscorr((P+1):end);
elseif strcmp(type, 'manual')   
    
    N_signal = length(rx_signal);
    Power = 1; %sqrt(preamble'*preamble).^2;
    
    if nargin < 4
        N_preamble = length(preamble);
        crosscorr = zeros(N_signal,1);
        for n = 1:(N_signal+N_preamble)        
            signal_causal = [zeros(N_preamble,1); rx_signal; zeros(N_preamble,1)];
            crosscorr(n) = sum( conj(preamble).*signal_causal(n:n+N_preamble-1) )/Power;            
        end
        crosscorr = crosscorr(N_preamble+1:end);
    
    elseif nargin == 4
        N_preamble = length(preamble)/2;    
        crosscorr = zeros(N_signal,1);
        %for n = 1:(N_signal+N_preamble)
        for n = 1:(N_signal-N+N_preamble)
            signal_causal = [zeros(N_preamble,1); rx_signal; zeros(N_preamble,1)];
            %crosscorr(n) = sum( conj(preamble).*signal_causal(n:n+N_preamble-1) )/Power;
            crosscorr(n) = sum( conj(preamble).*signal_causal([n:n+N_preamble-1 n+N:n+N+N_preamble-1]) )/Power;
        end
        crosscorr = crosscorr(N_preamble+1:end);
    elseif nargin == 5
        if noUW == 2
            N_preamble = length(preamble)/2;
            crosscorr = zeros(N_signal,1);
            %for n = 1:(N_signal+N_preamble)
            for n = 1:(N_signal-N+N_preamble)
                signal_causal = [zeros(N_preamble,1); rx_signal; zeros(N_preamble,1)];
                %crosscorr(n) = sum( conj(preamble).*signal_causal(n:n+N_preamble-1) )/Power;
                crosscorr(n) = sum( conj(preamble).*signal_causal([n:n+N_preamble-1 n+N:n+N+N_preamble-1]) )/Power;
            end
            crosscorr = crosscorr(N_preamble+1:end);
        elseif noUW == 3
            oneUW = preamble(1:end/2,:);            
            preamble = [preamble; oneUW];            
            N_preamble = length(preamble)/noUW;
            crosscorr = zeros(N_signal,1);
            for n = 1:(N_signal-2*N+N_preamble)
                signal_causal = [zeros(N_preamble,1); rx_signal; zeros(N_preamble,1)];
                crosscorr(n) = sum( conj(preamble).*signal_causal([n:n+N_preamble-1 n+N:n+N+N_preamble-1 n+2*N:n+2*N+N_preamble-1]) )/Power;
            end
            crosscorr = crosscorr(N_preamble+1:end);
        elseif noUW == 4            
            preamble = [preamble; preamble];            
            N_preamble = length(preamble)/noUW;
            crosscorr = zeros(N_signal,1);
            for n = 1:(N_signal-3*N+N_preamble)
                signal_causal = [zeros(N_preamble,1); rx_signal; zeros(N_preamble,1)];
                crosscorr(n) = sum( conj(preamble).*signal_causal([n:n+N_preamble-1 n+N:n+N+N_preamble-1 n+2*N:n+2*N+N_preamble-1 n+3*N:n+3*N+N_preamble-1]) )/Power;
            end
            crosscorr = crosscorr(N_preamble+1:end);
        elseif noUW == 8
            preamble = [preamble; preamble; preamble; preamble];
            N_preamble = length(preamble)/(noUW);
            crosscorr = zeros(N_signal,1);
            for n = 1:(N_signal-7*N+N_preamble)
                signal_causal = [zeros(N_preamble,1); rx_signal; zeros(N_preamble,1)];
                crosscorr(n) = sum( conj(preamble).*signal_causal([n:n+N_preamble-1 n+N:n+N+N_preamble-1 n+2*N:n+2*N+N_preamble-1 n+3*N:n+3*N+N_preamble-1 n+4*N:n+4*N+N_preamble-1 n+5*N:n+5*N+N_preamble-1 n+6*N:n+6*N+N_preamble-1 n+7*N:n+7*N+N_preamble-1]) )/Power;
            end
            crosscorr = crosscorr(N_preamble+1:end);
        end        
    end
end

end
