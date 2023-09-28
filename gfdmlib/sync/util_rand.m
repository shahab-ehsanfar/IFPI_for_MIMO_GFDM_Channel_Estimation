function randoms = util_rand(count, isComplex, seed)
% UTIL_RAND - Generate uniform random numbers, 1-dimensional
% use seed, if you want to start the sequence at a specific
% point. Otherwise, seed is assumed to be zero. Return uniform random
% values in [0, 1) in real and imaginary part.
% isComplex is a flag (i.e. 1 for complex, 0 for real (default))


% Implementation follows the linear congruential random generator.
% I know, the implmeentation is horrible, but, since the LVC random
% generator is not equal to the MATLAB one, even after seeding, we have to
% resort to our own random generator.

if nargin < 2
    isComplex = 0;
end

if nargin < 3
    seed = 0;
end



a = 1664525;
c = 1013904223;
state = seed * 500;
m = util_pow2(32);

randomsR = zeros(count, 1);
for i = 1:count
    state = mod(a * state + c, m);
    randomsR(i) = state / m;
end
randoms = randomsR;

if isComplex
    randomsI = zeros(count, 1);
    for i = 1:count
        state = mod(a * state + c, m);
        randomsI(i) = state / m;
    end
    randoms = randomsR + 1j*randomsI;
end
