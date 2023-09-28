function [Dout, Dp] = do_pilot_symbols_ofdm(p,Din,Delta_k,seq)
% Inserts pilot symbols into the first subsymbol of all subcarriers
% The pilot sequence generation is based on QPSK modulation accourding to
% 3GPP TS 36.211 Physical layer procedures, Release 8.9
%
% Din: D matrix after S/P
%
% Dout: Data with Pilot subsymbols
% Dp: Pilot subsymbols with zero on data positions
%
% Author: Shahab Ehsanfar

if nargin < 4
   seq = 'ZCH';
   if nargin < 3
       Delta_k = 1;
   end
end

if p.Mp == 2
    msetp = [1 p.Mp+p.Md];
else
    msetp = 1;
end

if strcmp(seq,'Rand')
    ZC1 = exp(2i*pi*randn(1,floor(p.K/Delta_k)+1)); %lteZadoffChuSeq(1, p.K+1);      %ones(1,p.K);%
else
    if mod(floor(p.K/Delta_k),2) == 0
        ZC1 = lteZadoffChuSeq(1, floor(p.K/Delta_k)*p.Mp+1);        
    else
        ZC1 = lteZadoffChuSeq(1, floor(p.K/Delta_k)*p.Mp+2);    
    end
end

Dd = Din; 
Dp = zeros(p.K,p.M);
for p_index = msetp
    if p_index == 1
        Din(p_index:Delta_k:end,1) = ZC1(1:p.Mp:p.Mp*floor(p.K/Delta_k));
    else
        Din(p_index:Delta_k:end,1) = ZC1(p.Mp:p.Mp:p.Mp*floor(p.K/Delta_k));
    end
    Dd(p_index:Delta_k:end,1) = 0;    
    Dp(p_index:Delta_k:end,1) = Din(p_index:Delta_k:end,1);
end

    Dout = Dd;
    



