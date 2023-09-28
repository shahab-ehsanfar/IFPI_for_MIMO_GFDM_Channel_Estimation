function [ rx_signal_lesscfo ] = ref_cfo_removal(preamble_starts, rx_signal, cfo_values, N )
% CFO Removal
%
% TU Dresden
% Shahab Ehsanfar




N_allblks = length(rx_signal);

cfo_values = [0; cfo_values];
preamble_starts = [preamble_starts preamble_starts(end)+N];

cfo_shifts = zeros(N_allblks,1);
for n = 0:N_allblks-1
    if n+1 < preamble_starts(1)+1
        blk = 0;
    elseif blk < length(preamble_starts)
        if n+1 > preamble_starts(blk+1)
            if blk+1 < length(preamble_starts)
                blk = blk + 1;
            end
        end    
    else
        error('CFO removal debugging needed')
    end    
%     x(n+1) = cfo_values(blk+1);
    cfo_shifts(n+1) = exp(-1j*2*pi*(n)*(cfo_values(blk+1))/N);
end


rx_signal_lesscfo = rx_signal.*cfo_shifts;


end

