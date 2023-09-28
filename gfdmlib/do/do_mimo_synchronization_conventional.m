function [ preamble_starts, cfo_values, rx_less_cfo, mean_cfo_hat ] = do_mimo_synchronization_conventional( config, rx_signal, mode )


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
   
   corr_stop = 100+2*Np+N;
   
%    corr = ref_sync_autocorrelate(distance,rx_signal, Ncp, Ncs);
%    
%    crosscorr_iT = zeros(length(rx_signal),nT);
%    for iT = 1:nT       
%        preamble = xp_iT(:,iT);       
%        crosscorr_iT(:,iT) = ref_sync_crosscorrelation(preamble,rx_signal,'manual');      
%    end   
%    crosscorr = prod(crosscorr_iT,2);

   
   for iR = 1:nR
    corr_iR(:,iR) = ref_sync_autocorrelate(Np,rx_signal(1:corr_stop,iR), Ncp, Ncs,'AC');
   end
   corr = sum(corr_iR,2); 
   
   cfo_flag = 1;
    if cfo_flag
        coarse_starts = ref_sync_peakfinding_relatedwork(distance, frame_length+N,corr);        
        
        cfo_values_iR = zeros(length(coarse_starts),nR);
        for iR = 1:nR
            cfo_values_iR(:,iR) = ref_frac_cfo_estimation_uw(coarse_starts, corr_iR(:,iR), N, Np,'primary');
        end
        cfo_values = mean(cfo_values_iR,2);
        
        for iR = 1:nR
            rx_less_cfo(:,iR) = ref_cfo_removal(coarse_starts, rx_signal(:,iR), cfo_values, N);
        end      
    else
       cfo_values = 0;
       rx_less_cfo = rx_signal;         
    end
        
   % CFO estimation & removal
%     cfo_hat = angle(corr(101))/(2*pi*nR*(-Np)/(N+Np));
    
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
   
   crosscorr_iT = zeros(corr_stop,nT*nR);
   for iR = 1:nR
       for iT = 1:nT
           preamble = xp_iT(:,iT);
           crosscorr_iT(:,(iR-1)*nT+iT) = ref_sync_crosscorrelation(preamble,rx_less_cfo(1:corr_stop,iR),'manual');
       end
   end
   crosscorr = sum(abs(crosscorr_iT),2);
     
   
   combined1 = crosscorr(1:length(corr)) / Np;
    
   
%    figure; plot(abs(combined1)./max(abs(combined1+1)),'ro');
%    hold on; plot(abs(corr)./max(abs(corr+1)),'m');
%    hold on; plot(abs(crosscorr)./max(abs(crosscorr))); legend('combined','autocorr','crosscorr')
   
   preamble_starts = ref_sync_peakfinding_relatedwork(distance, frame_length+N ,combined1);
   
   if length(preamble_starts) > 1
       preamble_starts = preamble_starts(1);
   end
   if length(cfo_values) > 1
       cfo_values = cfo_values(1);
   end
   
    
elseif strcmp(mode,'secondary')  
    xp_iT = config.preamble.xp_iT;
    xp_iT = [xp_iT(end/2+1:end,:); xp_iT(end/2+1:end,:)];
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
    corr_iR(:,iR) = ref_sync_autocorrelate(distance,rx_signal(:,iR), Ncp, Ncs,'AC','different_window',Np);
    end
    corr = sum(corr_iR,2);
%     corr = ref_sync_autocorrelate(distance,rx_signal, Ncp, Ncs,0,'manual',frame_length + config.preamble.N_preamble/2);
   
    cfo_flag = 1;
    if cfo_flag
        coarse_starts = ref_sync_peakfinding_relatedwork(distance, frame_length,corr);
        
        cfo_values_iR = zeros(length(coarse_starts),nR);
        for iR = 1:nR
            cfo_values_iR(:,iR) = ref_frac_cfo_estimation_uw(coarse_starts, corr_iR(:,iR), N, Np,'secondary');
        end
        cfo_values = mean(cfo_values_iR,2);
        
        mean_cfo_hat = repmat(mean(cfo_values),[length(coarse_starts) 1]);
        for iR = 1:nR
            rx_less_cfo(:,iR) = ref_cfo_removal(coarse_starts, rx_signal(:,iR), cfo_values, N);
        end        
    else
       cfo_values= 0; 
       rx_less_cfo = rx_signal;         
    end

%     CFO estimation & removal
%     cfo_hat = angle(corr(26))/(2*pi*nR*Np/(N+Np));
%     cfo_hat = angle(corr_iR(26,1))/(pi*(N+Np-1)/(Np));
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
        crosscorr_iT(:,(iR-1)*nT+iT) = ref_sync_crosscorrelation(xp_iT(:,iT),rx_less_cfo(:,iR),'manual',N);                 
    end
    end
    crosscorr = sum(abs(crosscorr_iT),2);
    
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
    
    combined = crosscorr(1:length(corr)) / distance;
    
%     figure; plot(abs(combined)./max(abs(combined+1)),'ro');
%     hold on; plot(abs(corr)./max(abs(corr)),'m');
%     hold on; plot(abs(crosscorr)./max(abs(crosscorr))); legend('combined','autocorr','crosscorr')
%     figure; plot(abs(combined)./max(abs(combined+1)));
    
    preamble_starts = ref_sync_peakfinding_relatedwork(distance, frame_length,combined);
%     [];
    if 0
        cfo_values_iR_new = zeros(length(preamble_starts),nR);
        for iR = 1:nR
            cfo_values_iR_new(:,iR) = ref_frac_cfo_estimation_uw(preamble_starts, corr_iR(:,iR), N, Np,'secondary');
        end
        updated_cfo_values = mean(cfo_values_iR_new,2);
        
        for iR = 1:nR
            rx_less_cfo(:,iR) = ref_cfo_removal(preamble_starts, rx_signal(:,iR), updated_cfo_values, N);
        end
        
        
        rx_signal = rx_less_cfo;
    end
    
    
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

