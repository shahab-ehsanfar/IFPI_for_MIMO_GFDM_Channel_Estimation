function [ preamble_starts ] = do_mimo_synchronization_relatedwork( config, rx_signal, mode )


if strcmp(mode,'primary')
   xp_iT = config.preamble.xp_iT;
   nT = size(xp_iT,2);
   nR = size(rx_signal,2);
    
   Ncp = 0;
   Ncs = 0;       
   
   frame_length = config.preamble.frame_length1;
   Np = config.Np;
   N = config.N;
   distance = Np;
%    corr = ref_sync_autocorrelate(distance,rx_signal, Ncp, Ncs);
%    
%    crosscorr_iT = zeros(length(rx_signal),nT);
%    for iT = 1:nT       
%        preamble = xp_iT(:,iT);       
%        crosscorr_iT(:,iT) = ref_sync_crosscorrelation(preamble,rx_signal,'manual');      
%    end   
%    crosscorr = prod(crosscorr_iT,2);

   
   for iR = 1:nR
    corr_iR(:,iR) = ref_sync_autocorrelate(Np,rx_signal(:,iR), Ncp, Ncs);
   end
   corr = sum(corr_iR,2);
   
       % CFO estimation & removal
%     cfo_hat = angle(corr(101))/(2*pi*nR*(-Np)/(N+Np));
%     
%     N_allblks = length(rx_signal);
%     cfo_shifts = zeros(N_allblks,1);
%     for n = 0:N_allblks-1
%         cfo_shifts(n+1) = exp(-1j*2*pi*(n)/(N+Np));
%     end    
%     rx_signal = (rx_signal).*(repmat(cfo_shifts,[1 nR]).^cfo_hat);

%     N_allblks = length(rx_signal);
% 
%     sto_shifts = zeros(N_allblks,1);
%     for n = 0:N_allblks-1
%         sto_shifts(n+1) = exp(+1j*2*pi*(n)/N_allblks);
%     end
%     sto = 0.2;
%     rx_signal = ifft(fft(rx_signal).*(repmat(sto_shifts,[1 nR]).^sto));
   
   crosscorr_iT = zeros(size(rx_signal,1),nT*nR);
   for iR = 1:nR
       for iT = 1:nT
           preamble = xp_iT(:,iT);
           crosscorr_iT(:,(iR-1)*nT+iT) = ref_sync_crosscorrelation(preamble,rx_signal(:,iR),'manual');
       end
   end
   crosscorr = sum(crosscorr_iT,2);
     
   
   combined1 = corr .* crosscorr(1:length(corr)) / Np;
    
   
%    figure; plot(abs(combined1)./max(abs(combined1+1)),'ro');
%    hold on; plot(abs(corr)./max(abs(corr+1)),'m');
%    hold on; plot(abs(crosscorr)./max(abs(crosscorr))); legend('combined','autocorr','crosscorr')
   
   preamble_starts = ref_sync_peakfinding_relatedwork(distance, frame_length ,combined1);
   
   if length(preamble_starts) > 1
       preamble_starts = preamble_starts(1);
   end
    
elseif strcmp(mode,'secondary')  
    xp_iT = config.preamble.xp_iT;
    nT = size(xp_iT,2);
    nR = size(rx_signal,2);
    noblk = config.noblk;
    
    Ncp = 0;
    Ncs = 0;
    
    
    N_data = config.payload.N_data;
    N = config.N;
    Np = config.Np;
    frame_length = config.preamble.frame_length2;
    distance = N;
          
    %corr = ref_sync_autocorrelate(distance,rx_signal, Ncp, Ncs,0,'different_window',frame_length + config.preamble.N_preamble/2);
    for iR = 1:nR
    corr_iR(:,iR) = ref_sync_autocorrelate(distance,rx_signal(:,iR), Ncp, Ncs,0,'different_window',Np);
    end
    corr = sum(corr_iR,2);
%     corr = ref_sync_autocorrelate(distance,rx_signal, Ncp, Ncs,0,'manual',frame_length + config.preamble.N_preamble/2);
    
%     % CFO estimation & removal
%     cfo_hat = angle(corr(33))/(2*pi*nR*Np/(N+Np));
%     
%     N_allblks = length(rx_signal);
%     cfo_shifts = zeros(N_allblks,1);
%     for n = 0:N_allblks-1
%         cfo_shifts(n+1) = exp(-1j*2*pi*(n)/(N+Np));
%     end    
%     rx_signal = (rx_signal).*(repmat(cfo_shifts,[1 nR]).^cfo_hat);

%    N_allblks = length(rx_signal);
% 
%     sto_shifts = zeros(N_allblks,1);
%     for n = 0:N_allblks-1
%         sto_shifts(n+1) = exp(+1j*2*pi*(n)/N_allblks);
%     end
%     sto = 0.2;
%     rx_signal = ifft(fft(rx_signal).*(repmat(sto_shifts,[1 nR]).^sto));

    crosscorr_iT = zeros(size(rx_signal,1)-Np,nT);
%     crosscorr_iT = zeros(length(rx_signal),nT); % Related Works
    for iR = 1:nR
    for iT = 1:nT
%         preamble = [xp_iT(end/2+1:end,iT); zeros(N_data,1); xp_iT(1:end/2,iT)];
        %preamble = [xp_iT(end/2+1:end,iT); xp_iT(1:end/2,iT)]; 
        crosscorr_iT(:,(iR-1)*nT+iT) = ref_sync_crosscorrelation(xp_iT(:,iT),rx_signal(:,iR),'manual',N);                 
    end
    end
    crosscorr = sum(crosscorr_iT,2);
    
%     figure; plot(abs(crosscorr_iT)./max(abs(crosscorr_iT(:))))
%     legend('$i_{\rm T} = 1$','$i_{\rm T} = 2$','$i_{\rm T} = 3$','$i_{\rm T} = 4$');
%     xlabel('$n$','interpreter','latex')
%     ylabel('$\Psi_{i_{\rm T}}[n]$','interpreter','latex')
%     title('SNR = 0 dB, $L = 9$','interpreter','latex')
%     
%     figure; plot(abs(prod(crosscorr_iT,2))./max(abs(prod(crosscorr_iT,2))))    
%     xlabel('$n$','interpreter','latex')
%     ylabel('$M_{i_{\rm R}, {\rm CC}}^{\rm (Pr)}[n]$','interpreter','latex')
%     title('SNR = 0 dB, $L = 9$','interpreter','latex')
%     
%     figure; plot(abs(sum(crosscorr_iT,2))./max(abs(sum(crosscorr_iT,2))))    
%     xlabel('$n$','interpreter','latex')
%     ylabel('$M_{i_{\rm R}, {\rm CC}}^{\rm (RW)}[n]$','interpreter','latex')
%     title('SNR = 0 dB, $L = 9$','interpreter','latex')
%     
%     for iT = 1:nT
%        xp_AC(:,iT) = xcorr(xp_iT(:,iT))/128;
%        figure(101); subplot(2,2,iT); plot(-127:127,abs(xp_AC(:,iT))); xlim([-127 127]); 
%     end    
%     figure; 
%     subplot(1,3,1)
%     plot(-127:127,abs(prod(xp_AC,2)))  
%     xlim([-127 127])
%     xlabel('$n$','interpreter','latex')
%     ylabel('$M_{i_{\rm R}, {\rm CC}}^{\rm (Pr)}[n]$','interpreter','latex')
%     subplot(1,3,2)
%     plot(-127:127,abs(sum(xp_AC,2)/4))    
%     xlim([-127 127])
%     xlabel('$n$','interpreter','latex')
%     ylabel('$M_{i_{\rm R}, {\rm CC}}^{\rm (RW)}[n]$','interpreter','latex')
%     subplot(1,3,3)
%     plot(-127:127,abs(prod(xp_AC,2).^4))  
%     xlim([-127 127])
%     xlabel('$n$','interpreter','latex')
%     ylabel('$M_{i_{\rm R}, {\rm CC}}^{\rm (Pr)}[n]$','interpreter','latex')
    
    combined = corr .* crosscorr(1:length(corr)) / distance;
    
%     figure; plot(abs(combined)./max(abs(combined+1)),'ro');
%     hold on; plot(abs(corr)./max(abs(corr)),'m');
%     hold on; plot(abs(crosscorr)./max(abs(crosscorr))); legend('combined','autocorr','crosscorr')
    
    
    preamble_starts = ref_sync_peakfinding_relatedwork(distance, frame_length,combined);
    
%     sync_counts = length(preamble_starts);
%     
%     while sync_counts < noblk        
%         blk_distances = zeros(sync_counts-1,1);
%         for i = 1:sync_counts-1
%            blk_distances(i) = preamble_starts(i+1) - preamble_starts(i);
%         end
%         missed_blk = find(blk_distances==max(blk_distances));
%         
%         correct_blk_distances = blk_distances;
%         correct_blk_distances(missed_blk) = [];
%         preamble_starts = [preamble_starts(1:missed_blk) preamble_starts(missed_blk)+floor(mean(correct_blk_distances)) preamble_starts(missed_blk+1:end)];
%         sync_counts = length(preamble_starts);
%     end
    
    
else
    error('Please select the synchronization mode (primary/secondary)')
    
end



end

