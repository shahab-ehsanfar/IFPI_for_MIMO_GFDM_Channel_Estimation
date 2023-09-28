function Dout = do_pilot_symbols_freq(p,Din)
% Inserts pilot symbols into the pilot subcarriers
% The pilot sequence generation is based on QPSK modulation accourding to
% 3GPP TS 36.211 Physical layer procedures, Release 8.9
%
% Din: D matrix after S/P
%
% Author: Shahab Ehsanfar

ksetp = (p.pf:p.pf:p.K); % pilot subcarriers

% Construct a PN object
h = commsrc.pn('Shift', 0);

% Set the output for p.M PN bits
set(h, 'NumBitsOut', p.M);

for p_index = ksetp
    
    pilots = generate(h); %round(rand(1,length(get_mset(p))));
    Din(p_index,:) = (1/sqrt(2))*(1-2*pilots) + 1j*(1/sqrt(2))*(1-2*circshift(pilots,1,2));
      
    assert(sum(pilots)>= 1);
    assert(sum(pilots)<= p.M-1);
    
end
Dout = Din;
