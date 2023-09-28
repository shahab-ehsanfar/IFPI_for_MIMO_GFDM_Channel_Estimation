% Copyright (c) 2014 TU Dresden
% All rights reserved.
% See accompanying license.txt for details.
%

function x = do_modulate(p, D, modType)
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
            if mod(p.M,2) == 0
                g = get_even_M_RC_filter(p.K, p.M, p.a);
            else
                g = get_transmitter_pulse(p);
            end
            g = g(round(1:K/L:end));
            G = fft(g);
        else
            G = p.cache.GL;
        end
        x = modulate_in_frequency(p, D, G);
    else
        if mod(p.M,2) == 0
            g = get_even_M_RC_filter(p.K, p.M, p.a);
        else
            g = get_transmitter_pulse(p);
        end
        x = modulate_in_time(p, D, g);
    end
else
    x = modulate_oqam(p, D);
end


end
