function [d_hat_soft, Heq] = calc_soft_equichan_v2(d_hat, config, R_dd_hat)
% Calculate the equivalent channel for soft demmaping and decoding



% R_dd_hat: a vector with diagonal elements of true R_dd_hat
          
Dd_indx = config.matrices.Dd_indx;



% GFDM
C = 1./(real(sqrt( (1./(R_dd_hat)) -1 )) );
d_hat(Dd_indx(:)) = d_hat(Dd_indx(:)).*C(Dd_indx(:));

Heq = C;%full(diag(equalizer*HA)); %diag(eye(size(equalizer*HA))); %


d_hat_soft = d_hat(Dd_indx(:));



Heq = Heq(Dd_indx(:)).';

end