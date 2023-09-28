% Copyright (c) 2014 TU Dresden
% All rights reserved.
% See accompanying license.txt for details.
%


function Dhat = demodulate_in_frequency_combine(p, G, x, pilot_set, norm_factor)
M = p.M; K = p.K;
L = length(G) / M;
LL = p.L;

indices = 1:M;
indices(pilot_set) = [];

pilot_subcar = p.Kset(1)+1:p.Delta_k:p.Kset(end)+1;
data_subcar = p.Kset+1;
data_subcar(pilot_subcar) = [];

Xhat = fft(x);
Dhat = zeros(K, M);
for k = pilot_subcar
    carrier = circshift(Xhat, ceil(L*M/2) - M*(k-1));
    carrier = fftshift(carrier(1:L*M));
    carrierMatched = carrier .* (G.*norm_factor);
    %dhat = downsample(ifft(carrierMatched), L);
    dhat = (sum(reshape(carrierMatched,M,L),2)/L);%M samples ifft
    dhat(indices) = ifft(dhat(indices));%*sqrt(M-1)/sqrt(M);
    Dhat(k,:) = dhat;
end
Dhat(pilot_subcar,pilot_set) = zeros(length(pilot_subcar),length(pilot_set)); %abs(Dhat(pilot_subcar,pilot_set)) > 0.01;

for k = data_subcar
    carrier = circshift(Xhat, ceil(L*M/2) - M*(k-1));
    carrier = fftshift(carrier(1:L*M));
    carrierMatched = carrier .* G;
    %dhat = downsample(ifft(carrierMatched), L);
    dhat = (sum(reshape(carrierMatched,M,L),2)/L);%M samples ifft
    dhat = ifft(dhat);%*sqrt(M-1)/sqrt(M);
    Dhat(k,:) = dhat;
end


end