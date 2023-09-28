function Dout = do_pilot_symbols_freq_a40(Np,iT,nT)
% Inserts pilot symbols into the pilot subcarriers
% The pilot sequence generation is based on QPSK modulation accourding to
% 3GPP TS 36.211 Physical layer procedures, Release 8.9
%
% Din: D matrix after S/P
%
% Author: Shahab Ehsanfar

N = Np*nT;
ksetp = 1:N; % pilot subcarriers

% Construct a PN object
h = commsrc.pn('Shift', 0);

% Set the output for p.M PN bits
set(h, 'NumBitsOut', 1);
pilots_bits = zeros(N,1);

for p_index = ksetp
    
    pilots = generate(h); %round(rand(1,length(get_mset(p))));
    %pilots_bits(p_index) = bi2de(str2double(sprintf('%d%d',pilots(1),pilots(2))));
    pilots_bits(p_index) = bi2de(pilots');
    %(1/sqrt(2))*(1-2*pilots) + 1j*(1/sqrt(2))*(1-2*circshift(pilots,1,2));
      
%     assert(sum(pilots)>= 1);
%     assert(sum(pilots)<= p.M-1);
    
end
Din = qammod(pilots_bits,2);
Dout = Din((iT-1)*Np+1:iT*Np);

