function [g_t, g_f] = get_even_M_RC_filter(K, M, alpha)
%K: number of samples per (sub)symbol, M: number of (sub)symbols
% alpha: roll-off factor
% Func: increasing function from -1 to 1;
% shift sampling shift
N = K*M;
%alpha;
if mod(M,2)== 0
    shift = 0.5;
else
    shift = 0;
end
f=(-1/2:1/N: (1/2-1/N))'+shift/N;
g_f = zeros(length(f),1);
ind = abs(f)<= (1-alpha)/(2*K);
g_f(ind) = 1;
ind = ((abs(f)> (1-alpha)/(2*K))& (abs(f) <= (1+alpha)/(2*K)));
f1 = f(ind);
g_f(ind)= 1/2*(1+sin(-pi*K/alpha*(abs(f1)-1/(2*K))));
g_f = sqrt(N)*g_f/norm(g_f);
g_f = fftshift(g_f);
g_t = ifft(g_f);
end
