function out = util_pow2(N)
% UTIL_POW2 - Calculate 2^n
out = 1;
for n = 0:N-1
    out = out * 2;
end
