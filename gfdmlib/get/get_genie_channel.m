function [ h_gains, h_gains_m, h_gains_all] = get_genie_channel(h,chan_typ,path_delays,numTaps,nT,nR,N,n_H)

if floor(path_delays) ~= path_delays % If path delays are fractional delays
    if strcmp(chan_typ,'time variant')
        h_gains = zeros(numTaps*nT,nR);
        h_gains_m = zeros(numTaps*nT,nR);
        for iR = 1:nR
            h_gains_iR = []; h_gains_iR_m = [];
            for iT = 1:nT
                h_i = h(iT,iR);
                
                h_genie = zeros(N,1); h_genie_m = zeros(N,1);
                for n = -N/2+1:N/2
                    for l = 1:length(path_delays)
                        h_genie(n+N/2) = h_genie(n+N/2) + h_i.PathGains(1,l)*sinc(path_delays(l) - n);
                        h_genie_m(n+N/2) = h_genie_m(n+N/2) + mean(h_i.PathGains(:,l))*sinc(path_delays(l) - n);
                    end
                end
                h_genie = circshift(h_genie,N/2+4); %%%% CHECK +4 in different scenarios
                h_genie(numTaps+1:end) = [];
                h_genie_m = circshift(h_genie_m,N/2+4); %%%% CHECK +4 in different scenarios
                h_genie_m(numTaps+1:end) = [];
                
                h_gains_iR = [h_gains_iR; h_genie];
                h_gains_iR_m = [h_gains_iR_m; h_genie_m];
            end
            
            h_gains(:,iR) = h_gains_iR;
            h_gains_m(:,iR) = h_gains_iR_m;
        end
    else
        h_gains = h;
        h_gains_m = h;
    end
else
    if strcmp(chan_typ,'time variant')
        h_gains{1} = zeros(numTaps*nT,nR); h_gains{2} = zeros(numTaps*nT,nR); h_gains{3} = zeros(numTaps*nT,nR);
        h_gains_m = zeros(numTaps*nT,nR);
        h_gains_all = zeros(numTaps,N,nT,nR);
        for iR = 1:nR
            h_gains_iR = [];
            for iT = 1:nT
                h_i = h(iT,iR);
                h_gains_iR = [h_gains_iR; h_i.PathGains.'];
                if size(h_i.PathGains,1) > N
                    h_gains_all(:,:,iT,iR) = (h_i.PathGains(end-N+1:end,:)).';
                elseif size(h_i.PathGains,1) == N
                    h_gains_all(:,:,iT,iR) = (h_i.PathGains(end-N+1:end,:)).';
                end
            end
            h_gains{1}(:,iR) = h_gains_iR(:,n_H(1)); % GFDM
            h_gains{2}(:,iR) = h_gains_iR(:,n_H(2)); % OFDM
            h_gains{3}(:,iR) = h_gains_iR(:,n_H(3)); % CP GFDM
            h_gains_m(:,iR) = mean(h_gains_iR,2);
        end
    else
        h_gains = h;
        h_gains_m = h;
    end
end

end

