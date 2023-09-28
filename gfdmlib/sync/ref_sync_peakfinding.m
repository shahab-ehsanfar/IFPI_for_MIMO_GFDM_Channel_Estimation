function preamble_starts = ref_sync_peakfinding(distance, frame_length, metric)

% Adaptive threshold 
no_of_thresholds = floor(length(metric)/frame_length);
T= zeros(no_of_thresholds,1);
for i= 1:no_of_thresholds
    T(i) = abs(ref_sync_estimateThreshold(metric(frame_length*(i-1)+1:frame_length*i),1));
end
const = ceil((length(metric)-no_of_thresholds*frame_length)/no_of_thresholds);
test = kron(T,ones(frame_length+const,1));
T = test(1:length(metric));

K = distance;
metric = abs(metric);
peak_pos = find(metric > T);
%if peak_pos(end) < length(metric) - K
    peak_pos = [peak_pos; length(metric)];
%end

preamble_starts = [];
p = 1;
while p < length(peak_pos)
    curr_pos = peak_pos(p);
    curr_peak = metric(peak_pos(p));

    while p <= length(peak_pos)
        p = p+1;
        if p > length(peak_pos); break; end
        next_pos = peak_pos(p);
        if ((next_pos - curr_pos) > K/2 || next_pos == peak_pos(end) )
            preamble_starts = [preamble_starts curr_pos];
            p = find(peak_pos >= curr_pos + frame_length-K/2, 1);
            if isempty(p); p = length(peak_pos); end;
            break;
        end
        next_peak = metric(next_pos);
        if (next_peak > curr_peak); break; end;
    end
end
