function [ bitErr_gen, ser, mi ] = do_detection_test(p,N,nT,b,bc,dd,Y,Heq,snr,PCCC,code_interleaver, bicm_interleaver,Ainv_F_mimo,F_mimo_A,p2)


if nargin < 15 || p2.M ~= 1
    EsNo = 10^(snr/10);
    [Xhat,Q, QHeq] = do_chan_equalization(N,Heq,Y,snr,nT);
    Dhat = zeros(p.K,p.M,nT);
    [C_genie, E_h] = calc_softdec_equichan(p,Ainv_F_mimo,F_mimo_A,QHeq,Q,nT);
else
    [Xhat,Q, QHeq] = do_chan_equalization(N,Heq,Y,snr,nT,p2);
    Dhat = zeros(p2.K,1,nT);
    [C_genie, E_h] = calc_softdec_equichan(p2,Ainv_F_mimo,F_mimo_A,QHeq,Q,nT);
    EsNo = 1;% 10^(snr/10);
end



H_equiv = C_genie.*E_h;
N_d = length(H_equiv)/nT;

Delta_k = p.Delta_k;


dhat = zeros(p.Kon*p.M-p.Mp*(p.Kon/Delta_k),nT);

bit_llr = zeros(size(bicm_interleaver,2),nT);
if (isfield(PCCC,'hard') && PCCC.hard == true)
    bh = bit_llr; 
else
    bh = zeros(size(code_interleaver,2),nT); 
end

for iT = 1:nT
   
   if nargin < 15
       Dhat(:,:,iT) = do_demodulate_combine(p, ifft_u(Xhat(N*(iT-1)+1:iT*N)),'ZF',1:nT);
       dhat(:,iT) = do_unmap_p(p, Dhat(:,:,iT));
   elseif p2.M == 1
       Dhat(:,:,iT) = do_demodulate(p2, ifft_u(Xhat(N*(iT-1)+1:iT*N)),'ZF');
       Dhat_ofdm_var = reshape(Dhat(:,:,iT),p.M,p.K).';
       dhat(:,iT) = do_unmap_p(p,Dhat_ofdm_var);      
%        d_set = 1:N;
%        if p.Mp == 1
%            d_set(1:p.Delta_k*p.M:N) = [];
%        elseif p.Mp == 2
%            d_set([1:p.Delta_k*p.M:N 2:p.Delta_k*p.M:N]) = [];
%        end
%        dhat(:,iT) = Dhat(d_set,:,iT);
       
   else
       Dhat(:,:,iT) = do_demodulate(p2, ifft_u(Xhat(N*(iT-1)+1:iT*N)),'ZF');
       dhat(:,iT) = do_unmap_p(p,Dhat(:,:,iT));
   end  
   E_h = ones(size(dhat(:,iT))).';
    n0 = (iT-1)*N_d;
    EsNo = 10^(snr/10);
   [bh(:,iT), bit_llr(:,iT)] = do_decode_pccc(dhat(:,iT), b(:,iT).', EsNo, E_h, PCCC, code_interleaver(iT,:), bicm_interleaver(iT,:));
   %[bh(:,iT), bit_llr(:,iT)] = do_decode_pccc((C_genie(n0+1:n0+N_d).').*dhat(:,iT), b(:,iT).', EsNo, H_equiv(n0+1:n0+N_d), PCCC, code_interleaver(iT,:), bicm_interleaver(iT,:));
   %[bh(:,iT), bit_llr(:,iT)] = do_decode_pccc_loopedDemod((C_genie(n0+1:n0+N_d).').*dhat(:,iT), b(:,iT).', EsNo, H_equiv(n0+1:n0+N_d), PCCC, code_interleaver(iT,:), bicm_interleaver(iT,:));
end

mi = mutualinfo(bit_llr(:),bc(:));

if (isfield(PCCC,'hard') && PCCC.hard == true)    
    bitErr_gen = sum( bh(:) ~= bc(:) );
else
    bitErr_gen = sum( bh(:) ~= b(:) );
end

s = do_qamdemodulate(dd(:),p.mu);
sh = do_qamdemodulate(dhat(:)./(E_h.'),p.mu);

ser = sum( sh ~= s );

% sh_uneq = do_qamdemodulate(dhat(:),p.mu);
% bitErr_gen = sum( sh_uneq ~= s );

end

