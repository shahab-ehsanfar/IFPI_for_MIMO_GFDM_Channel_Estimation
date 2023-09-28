% Copyright (c) 2014 TU Dresden
% All rights reserved.
% See accompanying license.txt for details.
%


function X = fft_u(x)
% Unitary DFT
%
% Normalizes with respect to sqrt(N)


X = fft(x)/sqrt(size(x,1));


