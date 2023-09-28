function [d_hat_soft, Heq] = do_detect_soft_sym(i_blk, p, conf, pilots_conf, Y, snr, vH_hat, R_hh_err, N, nT, nR, ofdm_flag, d_hat_soft, Heq)

pilots_conf.blk = i_blk;
 if ofdm_flag == true       
     [d_hat_ofdm{i_blk}, ~, equalizer_ofdm, R_dd_hat_ofdm, HA_ofdm] = do_CEPIC_equalization( conf, pilots_conf, snr, reshape(vH_hat, [N nT*nR]), R_hh_err, Y(:) );
     [d_hat_soft{i_blk}, Heq{i_blk}] = calc_soft_equichan_LMMSE(d_hat_ofdm{i_blk}, pilots_conf, equalizer_ofdm, R_dd_hat_ofdm, HA_ofdm, nT, ofdm_flag, [p.K p.M]);          
 else
     [d_hat{i_blk}, ~, equalizer, R_dd_hat, HA] = do_CEPIC_equalization( conf, pilots_conf, snr, reshape(vH_hat, [N nT*nR]), R_hh_err, Y(:,i_blk));     
     [d_hat_soft{i_blk}, Heq{i_blk}] = calc_soft_equichan_LMMSE(d_hat{i_blk}, pilots_conf, equalizer, R_dd_hat, HA, nT);
 end
 

end