% Copyright (c) 2014 TU Dresden
% All rights reserved.
% See accompanying license.txt for details.
%


function [rise,fall]=get_ramp(M,a,type)
m=(0:(M-1));
R=((0:(M-1))'-M/2-eps)/(a*M)+1/2;R(R<0)=0;R(R>1)=1;F=1-R;% ramp rise
R4th=R.^4.*(35 - 84*R+70*R.^2-20*R.^3);F4th=1-R4th;F4th=1-R4th;% 4th order rise
%falling ramp
switch type
    case 0, rise=R.^.5;                      fall=(1-rise.^2).^.5;%root ramp 
    case 1, rise=R;                          fall=1-rise;         %ramp 
    case 2, rise=R4th.^.5;                   fall=(1-rise.^2).^.5;%root 4th order
    case 3, rise=R4th;                       fall=1-rise;         %4th order
    case 4, rise=(1/2*(cos(F*pi)+1)).^.5;    fall=(1-rise.^2).^.5;%root raised cossine
    case 5, rise=1/2*(cos(F*pi)+1);          fall=1-rise;         %raised cossine
    case 6, rise=(1/2*(cos(F4th*pi)+1)).^.5; fall=(1-rise.^2).^.5;%rrc 4th oder
    case 7, rise=1/2*(cos(F4th*pi)+1);       fall=1-rise;         %rc 4th order
    case 8, rise=1/2*(exp(-1i*F*pi)+1);      fall=1-rise;         %xia
    case 9, rise=1/2*(exp(-1i*F4th*pi)+1);   fall=1-rise;         %xia 4th order
    case 10, % Dolph-Chebyshev
        N = 2*M;
        alpha = 6;
        Beta = cosh((1/N)*acosh(10^alpha));
        n = (0:(N-1));
        W = cos(N*acos(Beta*cos(pi*n./N)))./(cosh(N*acosh(Beta)));
        w = ifft(W).';  
        rise = w(end/2+1:end)/max(w); fall = w(1:end/2)/max(w); 
    case 11, % Blackman
        alpha = 0.16;
        N = 2*M;
        n = (0:(N-1));
        a0 = (1-alpha)/2; a1 = 0.5; a2 = alpha/2;
        w = (a0 - a1*cos(2*pi*n./(N-1))+a2*cos(4*pi*n./(N-1))).';
        fall = w(end/2+1:end)/max(w); rise = w(1:end/2)/max(w); 
end
