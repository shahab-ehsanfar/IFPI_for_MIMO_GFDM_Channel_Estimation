function [on_bins, on_bins1, off_bins, on_bins_single] = get_on_bins_ofdm(p, nT,nR)

N= p.K;

on_bins_single = zeros(N,1);

on_bins_single(p.Kset+1) = 1; 

on_bins_single(p.Kpset+1) = 0; % Kpset: Kset of pilots

p_ofdm.Kset = find(on_bins_single) - 1;
p_ofdm.Kon = length(p_ofdm.Kset);

on_bins_single = p.Kset+1;


on_bins1 = [];
for iT = 0:nT-1
   on_bins1 = [on_bins1; (on_bins_single+iT*N)]; 
end

on_bins = [];
for iR =  0:nR-1
    on_bins = [on_bins; (on_bins1+iR*nT*N)];
end



% Xp_iT(abs(Xp_iT) < 1e-9 ) = 0;
% bin_start = find(abs(Xp_iT(:,1)),1); 
% bin_end = find(abs(Xp_iT(:,1)),1,'last'); 
% noChannels = nT*nR; 
% on_bins = [];
% for i_Ch = 0:noChannels-1
%     on_bins = [on_bins i_Ch*N+bin_start:i_Ch*N+bin_end];
% end
% 
% %if Koff > 10
%     on_bins1 = [];
%     for i_Ch = 0:nR-1
%         on_bins1 = [on_bins1 i_Ch*N+bin_start:i_Ch*N+bin_end]; % on_bins1 = [on_bins1 i_Ch*N+bin_start-p.M:i_Ch*N+bin_end+p.M];
%         if i_Ch == 0
%             on_bins_single = on_bins1;
%         end
%     end
% %else
% %     on_bins1 = [];
% %     for i_Ch = 0:nR-1
% %         on_bins1 = [on_bins1 i_Ch*N+bin_start:i_Ch*N+bin_end];
% %         if i_Ch == 0
% %             on_bins_single = on_bins1;
% %         end
% %     end
% %end
off_bins = 1:N*nR;
off_bins(on_bins1) = [];


end
