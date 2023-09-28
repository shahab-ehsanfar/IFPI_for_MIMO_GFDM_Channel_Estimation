function [d_hat_soft, Heq] = calc_soft_equichan_LMMSE(d_hat, pilots_conf, equalizer, R_dd_hat, HA, nT, ofdm_flag, KM)
% Calculate the equivalent channel for soft demmaping and decoding

% Component-wise unbiased LMMSE
if pilots_conf.blk == 1
    Dd_indx = logical(pilots_conf.Dd1_indx);
else
    Dd_indx = logical(pilots_conf.Dd_indx);
end


if nargin > 6
    if ofdm_flag == true  % OFDM
                
        Es = 1;
        C = 1./(real(sqrt( (1./diag(R_dd_hat))  -1 )));  % Remove the Es = mean(dd)
        d_hat_full = zeros(prod(KM)*nT,1);
        d_hat_full(Dd_indx(:)) = d_hat(Dd_indx(:)).*C(Dd_indx(:));
        
        Dd_hat_ofdm = zeros([KM nT]);
        C_ofdm = zeros([KM nT]);
        Dd_indx_ofdm = zeros([KM nT]);
        d_hat_nT = reshape(d_hat_full, [prod(KM) nT]);
        C_nT = reshape(C, [prod(KM) nT]);
        
        Dd_hat = zeros(prod(KM),nT); C_demapped = zeros(prod(KM),nT);
        for iT = 1:nT
            Dd_hat_ofdm(:,:,iT) = reshape(d_hat_nT(:,iT), [KM(2) KM(1)] ).';            
            C_ofdm(:,:,iT) = reshape(C_nT(:,iT), [KM(2) KM(1)] ).';
            Dd_indx_ofdm(:,:,iT) = reshape(Dd_indx(:,:,iT), [KM(2) KM(1)]).';
            
            Dd_hat(:,iT) = reshape( Dd_hat_ofdm(:,:,iT) , [prod(KM) 1]);
            C_demapped(:,iT) = reshape( C_ofdm(:,:,iT) , [prod(KM) 1]);
            Dd_indx(:,:,iT) = reshape( Dd_indx_ofdm(:,:,iT) , [prod(KM) 1]);
        end
        d_hat = Dd_hat(:);
        C_final = C_demapped(:);
        
        Heq = C_final;
        
        
        
    end
    
else  % GFDM
    C = 1./(real(sqrt( (1./diag(R_dd_hat)) -1 )) );
    d_hat(Dd_indx(:)) = d_hat(Dd_indx(:)).*C(Dd_indx(:));
    
    Heq = C;%full(diag(equalizer*HA)); %diag(eye(size(equalizer*HA))); %
end

d_hat_soft = d_hat(Dd_indx(:));



Heq = Heq(Dd_indx(:)).';

end