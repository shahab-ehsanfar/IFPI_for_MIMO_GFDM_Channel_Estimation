% Copyright (c) 2014 TU Dresden
% All rights reserved.
% See accompanying license.txt for details.
%

function x = do_modulate_precode(p, D, modType, pilot_set)
% Modulate the GFDM data matrix with the GFDM system
%
% Uses efficient FFT implementation.
% \param D KxM-matrix containing the data to be transmitted at each
%          subcarrier and symbol.
% \param modType is optional and defines
%    - F: modulate in frequency 
%    - T: modulate in time (default)

L = p.L; K = p.K; M = p.M; oQAM = p.oQAM;
N = M*K;


if ~oQAM
    if nargin > 2 && modType=='F'
        if ~isfield(p, 'cache.GL')
            g = get_transmitter_pulse(p);
            g_dwnsmpl = g(round(1:K/L:end));
            G = fft(g_dwnsmpl);
            G(:,2) = G;
        else
            G = p.cache.GL;
        end
        G_precode = fft(g); G_precode(1:2) = [0;0];
        g_precode = ifft(G_precode);
        norm_factor = sqrt(sum(abs(g_precode).^2));
        G(:,2) = G(:,2)/norm_factor; % Renormalize filter at pilot subcarriers    
        x = modulate_in_frequency_precode(p, D, G, pilot_set);
    else
        g = get_transmitter_pulse(p);
        x = modulate_in_time(p, D, g);
    end
else
    x = modulate_oqam(p, D);
end


end
