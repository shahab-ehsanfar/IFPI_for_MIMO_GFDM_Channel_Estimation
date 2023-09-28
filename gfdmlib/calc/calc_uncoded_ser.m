function [ser] = calc_uncoded_ser(symbols1, symbols2, d_hat, Heq, p)

symb_hat1 = do_qamdemodulate(d_hat{1}./Heq{1}.' ,p.mu);
symb_hat2 = do_qamdemodulate(d_hat{2}./Heq{2}.' ,p.mu);

ser = mean( [symbols1; symbols2] ~= [symb_hat1; symb_hat2] ); 