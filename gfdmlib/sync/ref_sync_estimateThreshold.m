function T = ref_sync_estimateThreshold(metric,constant)
%constant = 0.5;%0.05;
T = constant*mean(metric) + constant*std(metric);
