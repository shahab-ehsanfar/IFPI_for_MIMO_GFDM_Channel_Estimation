function [ cfo_values ] = ref_frac_cfo_estimation_uw(preamble_starts, autocorrelation, N, Np, mode)
% CFO Estimation for Unique Word transmission
%
% TU Dresden
% Shahab Ehsanfar

if strcmp(mode,'primary') % Integer CFO upto 0.5*(N/Np)
    CFO_fr_Np = angle(autocorrelation(preamble_starts))/(2*pi);
    cfo_values = CFO_fr_Np*(N/(Np)); 
    
elseif strcmp(mode,'secondary') % Only Fractional CFO upto 0.5 subcarrier spacing
    no_detected_blks = length(preamble_starts);
    cfo_values = zeros(no_detected_blks,1);
    for i = 1:no_detected_blks
        cfo_values(i) = angle(autocorrelation(preamble_starts(i)))/(2*pi);
    end
elseif strcmp(mode,'secondaryCirc') % Only Fractional CFO upto 0.25 subcarrier spacing
    no_detected_blks = length(preamble_starts);
    cfo_values = zeros(no_detected_blks,1);
    for i = 1:no_detected_blks
        cfo_values(i) = angle(autocorrelation(preamble_starts(i)))/(4*pi);
    end    
end

end

