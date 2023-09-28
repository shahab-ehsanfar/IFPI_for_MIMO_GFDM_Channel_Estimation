% Copyright (c) 2014 TU Dresden
% All rights reserved.
% See accompanying license.txt for details.
%


function A = get_transmitter_matrix_precode(p)


M = p.M; K = p.K; L = p.L;

g = get_transmitter_pulse(p);
g_dwnsmpl = g(round(1:K/L:end));
G = fft(g_dwnsmpl);

G_precode = fft(g); G_precode(1:2) = [0;0];
g_precode = ifft(G_precode);
norm_factor = sqrt(sum(abs(g_precode).^2));
G = G/norm_factor;

dftsize_Md.M = M - p.Mp; dftsize_Md.K = 1;
P_0 = [eye(M) zeros(M,M*(K-1)); zeros(M,M*(K-1)) eye(M)].';

Gamma = diag(G);

R = [eye(M) eye(M)].';

F_Md = get_DFT_matrix(dftsize_Md);

pilot_scale = sqrt(M)./(max(L*G)*sqrt(M-1));

Af_new = [];
for k = 0:K-1
    W_k = blkdiag(pilot_scale,pilot_scale,F_Md*sqrt(p.Md));
        
    P_k = circshift(P_0,M*k);
    
    A_k = P_k*Gamma*R*W_k;
    
    Af_new = [Af_new A_k];
end

F = get_DFT_matrix_nounit(p);

A = (K/L)*F'*Af_new/(M*K);


% for D transposed version
A_transposed = zeros(M*K);
for n = 1:M*K
    for m = 1:M
        A_transposed(n,1+(m-1)*K:m*K) = A(n,m:M:M*K);
    end
end

clear A
A = A_transposed;

end
