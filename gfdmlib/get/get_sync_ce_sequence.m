function [ xp_iT ] = get_sync_ce_sequence( N_preamble, iT, type, PrimeNr, no_skip_samples)
% N_preamble: length of (double sequence) preamble
% PrimeNr: The Prime number for length of Polyphase sequence P^2


if strcmp(type,'upchirp')
    N_half = N_preamble/4;
    N_chirps = 1;
    xp0_half = exp(1i.*pi.*N_chirps/N_half.*(0:N_half-1).^2).'./sqrt(N_half); % UP CHIRP
    %             xp0 = 5*[xp0_half; conj(xp0_half)];
    
    %             N_quart = N/4;
    %             N_chirps = 1;
    %             xp0_quart = exp(1i.*pi.*N_chirps/N_quart.*(0:N_quart-1).^2).'./sqrt(N_quart); % UP CHIRP
    
    xp0_dbl = ifft(upsample(fft(xp0_half),2));
    
    %             xp0 = 5*[xp0_half; (xp0_half)];
    %xp0 = 5*[xp0_dbl(1:end/2); conj(xp0_dbl(1:end/2))];
    xp0 = 5*xp0_dbl;
    
    %             numPaths = 4;
    %             N_half = 16;
    %             if mod(numPaths,2) == 0
    %                 nu = 2;
    %             else
    %                 nu = 1;
    %             end
    %
    %             %if (iT == 1 || iT ==3)
    %                 xp1 = zeros(N_half+numPaths-1,1);
    %                 for n = -numPaths+1:(N_half-1)
    %                     xp1(n+numPaths) = exp(1i*2*pi*(1-1)*n/(nT*numPaths) )*exp(1i*2*pi*n*(n+nu)/numPaths);
    %                 end
    %                 xp1_long = 5*[xp1; (xp1)];
    %             %elseif (iT == 2 || iT ==4)
    %                 xp2 = zeros(N_half+numPaths-1,1);
    %                 for n = -numPaths+1:(N_half-1)
    %                     xp2(n+numPaths) = exp(1i*2*pi*(2-1)*n/(nT*numPaths) )*exp(1i*2*pi*n*(n+nu)/numPaths);
    %                 end
    %                 xp2_long = 5*[xp2; (xp2)];
    %             %end
    %
    %             % orthogonality condition
    %             sum = zeros(numPaths,numPaths);
    %             for k = 0:numPaths-1
    %                 for l = 0:numPaths-1
    %                     for n = 0:N_half-1
    %                         product = (xp1(n+numPaths)*conj(xp1(n+numPaths)))*exp(1i*2*pi*(k-l)*n/N);
    %                         sum(k+1,l+1) = sum(k+1,l+1) + product/N_half;
    %                     end
    %                 end
    %             end
    %             figure; mesh(abs(sum))
    
    if (iT == 1 || iT == 3)
        % xp_iT(:,iT) =  [xp0(1:end/2); xp0]; %(fft(xp0).*(shifts.^(l_0*(iT-1))));
        xp_iT =  [xp0; xp0];
        %xp_iT(:,iT) =  circshift(do_modulate(p, Dp(:,:,iT)),l_0*(iT-1));
        %else
        %xp_iT(:,iT) =  zeros(size(do_modulate(p, Dp(:,:,iT)),l_0*(iT-1)-1));
        
        %xp_iT(:,iT) = xp1_long;
    else
        xp2 = 5*ifft(circshift(fft(xp0_dbl),1));
        % xp_iT(:,iT) =  ([xp2(end/2+1:end); xp2]); %[xp2; conj(xp2)];
        xp_iT = ([xp2; xp2]);
        %xp_iT(:,iT) =  conj(flipud(xp0)); %ifft(fft(xp0).*(shifts.^(l_0*(iT-1))));
        
        %                  xp_iT(:,iT) =  circshift(conj(flipud(xp0)),1);
        %xp_iT(:,iT) =  ifft(fft(do_modulate(p, Dp(:,:,1))).*(shifts.^(l_0*(iT-1))));
        %             elseif iT == 3
        %                 xp_iT(:,iT) =  ifft(fft(xp0).*(shifts.^(l_0*(iT-1))));
        %                 %xp_iT(:,iT) =  ifft(fft(do_modulate(p, Dp(:,:,1))).*(shifts.^(l_0*(iT-1))));
        %             elseif iT == 4
        %                 xp_iT(:,iT) =  ifft(fft(xp0).*(shifts.^(l_0*(iT-1))));
        %                 %xp_iT(:,iT) =  ifft(fft(do_modulate(p, Dp(:,:,1))).*(shifts.^(l_0*(iT-1))));
        %xp_iT(:,iT) = xp2_long;
    end
elseif strcmp(type,'upchirp4')
    nT =4;
    N_half = N_preamble/(2*nT);
    N_chirps = 1;
    xp0_half = exp(1i.*pi.*N_chirps/N_half.*(0:N_half-1).^2).'./sqrt(N_half); % UP CHIRP   
    
    xp0_dbl = ifft(upsample(fft(xp0_half),4));
    
    xp0 = 5*xp0_dbl;
    
    
    
    if iT == 1         
        xp1 = 5*ifft(circshift(fft(xp0_dbl),0));        
        xp_iT =  [xp1; xp1];        
    elseif iT == 2       
        xp2 = 5*ifft(circshift(fft(xp0_dbl),1));        
        xp_iT = ([xp2; xp2]);
    elseif iT == 3
        xp3 = 5*ifft(circshift(fft(xp0_dbl),2));        
        xp_iT = ([xp3; xp3]);
    elseif iT == 4
        xp4 = 5*ifft(circshift(fft(xp0_dbl),3));        
        xp_iT = ([xp4; xp4]);
    end
elseif strcmp(type,'polyphase')
    
    
    P = PrimeNr; %sqrt(N_preamble/2); 
%     assert(isprime(P),'Wrong input: P must be a prime number');
    
    % Works well with N = 128, P^2 = 144, distance 126
    %N_preamble = 144;
    xp_iT_polyph = zeros(P^2,1);
    for n0 = 0:P-1
        for n1 = 0:P-1
            n = n0*P+n1;
            xp_iT_polyph(n+1) = exp(1i*2*pi*iT*n0*n1/P);
        end        
    end
    
    xp_iT_single = ifft_u([zeros(floor((N_preamble/2 - P^2)/2),1); fft_u(xp_iT_polyph); zeros(ceil((N_preamble/2 - P^2)/2),1)]);
    
    xp_iT_single(1:no_skip_samples) = [];
    
%     xp_iT = [conj(xp_iT_single); xp_iT_single];
    xp_iT = [(xp_iT_single); xp_iT_single];
    
%     xp_iT(end/2+1:end) = xp_iT(1:end/2);
    
      
elseif strcmp(type, 'randperm4')    
    nT =4;
    N_half = N_preamble;%/(2*nT);
    N_chirps = 1;
    xp0_half = exp(1i.*pi.*N_chirps/N_half.*(0:N_half-1).^2).'./sqrt(N_half); % UP CHIRP   
    
    xp0_dbl = ifft(upsample(fft(xp0_half),4));
    
    xp0 = 5*xp0_dbl;
    
    xp_bins = [];
    for iT = 1:N_preamble/nT
       xp_bins = [xp_bins; (iT-1)*nT+randperm(nT)']; 
    end
    x_test = fft(xp0_half);
    x_test(xp_bins(2:4:end)) = 0;
    x_test(xp_bins(3:4:end)) = 0;
    x_test(xp_bins(4:4:end)) = 0;
    figure; plot(abs(x_test)); figure; plot(abs(ifft(x_test)))
    xp0_new = ifft(x_test);
    
    % Manual cross correlation (convolution)
    signal = [zeros(200,1); xp0_new; zeros(100,1); xp0_new];
    crosscorr = zeros(length(signal),1);
    for n = 1:(length(signal)+N_preamble+1)
       signal_causal = [zeros(N_preamble,1); signal; zeros(N_preamble,1)];       
       crosscorr(n) = sum( conj(xp0_new).*signal_causal(n:n+N_preamble-1) );       
    end
    crosscorr = crosscorr/max(crosscorr);
    figure; plot(abs(crosscorr(N_preamble+1:end)))
    
    %%%%%%%%%%
%     TEST FOR CROSS CORRELATION OVER MULTI ANTENNAS

%     crosscorr_iT = zeros(length(signal),nT);
%     for iT = 1:nT    
%     xp0_new= xp_iT(:,iT);
%     signal = [zeros(200,1); xp0_new; zeros(100,1); xp0_new];
%     crosscorr = zeros(length(signal),1);
%     for n = 1:(length(signal)+N_preamble)
%        signal_causal = [zeros(N_preamble,1); signal; zeros(N_preamble,1)];       
%        crosscorr(n) = sum( conj(xp0_new).*signal_causal(n:n+N_preamble-1) );       
%     end
%     crosscorr_iT(:,iT) = crosscorr(N_preamble+1:end);
%     end
%     
%     figure; plot(abs(prod(crosscorr_iT,2)))
%     
end


end

