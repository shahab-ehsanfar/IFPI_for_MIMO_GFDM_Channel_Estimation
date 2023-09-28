% Copyright (c) 2014 TU Dresden
% All rights reserved.
% See accompanying license.txt for details.
%


function x = modulate_in_frequency_precode(p, D, G, pilot_set)
K = p.K; M = p.M;
N = M*K;
L = length(G(:,1)) / M;

% pp = p; pp.M = p.Md;
% g = get_transmitter_pulse(pp);
% g = g(round(1:K/L:end));
% G = fft(g);
% G = [1; 1; G(1:end/2); 0; 0; G(end/2+1:end)];

block0 = D.';


indices = 1:size(block0,1);
indices(pilot_set) = [];

pilot_subcar = p.Kset(1)+1:p.Delta_k:p.Kset(end)+1;
data_subcar = p.Kset+1;
%data_subcar(pilot_subcar) = [];
data_subcar(pilot_subcar-p.Kset(1)) = [];

% FFT over data only
block = block0;
%block(indices,:) = fft(block0(indices,:), [], 1);
block(indices,pilot_subcar) = fft(block(indices,pilot_subcar), [], 1);
block(:,data_subcar) = fft(block(:,data_subcar), [], 1);

% The following line makes sure pilots signal has energy equal to 1
%block(pilot_set,pilot_subcar) = block(pilot_set,pilot_subcar)./max(L*G);
block(pilot_set,pilot_subcar) = block(pilot_set,pilot_subcar)*( sqrt(M) / (sqrt(M-1)*L*max(G(:,2))) );

DD = repmat(block, L, 1);
X = zeros(N,1);

for k=1:K
    carrier = zeros(N,1);
    if any(abs(k - pilot_subcar)<1e-10)
        carrier(1:L*M) = fftshift(DD(:,k) .* G(:,2));
    else
        carrier(1:L*M) = fftshift(DD(:,k) .* G(:,1));
    end
    carrier = circshift(carrier, -floor(L*M/2) + M*(k-1));
%    hold on; plot(abs((0.58*K/L)*carrier),'r')
    X = X + carrier;
end

X = (K/L) * X;


x = ifft(X);


end
