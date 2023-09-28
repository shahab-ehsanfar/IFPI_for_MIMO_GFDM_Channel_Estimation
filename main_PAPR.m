% Copyright (c) 2014 TU Dresden
% All rights reserved.
% See accompanying license.txt for details.
%
clear variables

snr = -5:5:35;
p=get_defaultGFDM('OFDM');
p.pulse='rc_fd';
p.a=0.1;%roll-off factor for the RRC and RC filters.
p.K=4;%Number of useful subcarriers.
p.L=2;%Number of overlapping subcarrier for frequency domain GFDM implementation..
p.Kset = 0:p.K(end)-1;
p.Kon = length(p.Kset);
%p.Kon=floor(p.K/2)+1;%Number of actives subcarriers.
p.Mp = 2; % Number of pilot subsymbols
p.Md = 19; % Number of data subsymbols 
p.M= p.Md + p.Mp; % Number of subsymbols.
p.mu=4; % 16-QAM.
p.Ncp = p.K; % Number of CP samples
p.b = 0; %p.K/4; % Number of window raise/fall samples
p.window = 'rc';
%p.window = 'rc';
p.mode = 'precode';

its = 20000;

% Currently pilot insertion is adapted for only up to 2x2 MIMO channel
nT = 1; % Number of Tx antennas (up to 2)
nR = 1; % Number of Rx antennas (up to 2)

p.Delta_k = 2; % pilot subcarrier spacing
numPaths = 24; % number of non-zero channel taps

p_ofdm = p;
p_ofdm.M = 1;
p_ofdm.K = p.M*p.K;
p_ofdm.a = 0;
p_ofdm.L = p_ofdm.K;
p_ofdm.pulse = 'rc_td';

N = p.M*p.K;
Delta_k = p.Delta_k;


errorsZF_gen = zeros(length(snr), its);
errorsZF_ls = zeros(length(snr), its);
errorsZF_lmmse = zeros(length(snr), its);


P = get_power_delay_profile(numPaths,nT,nR); 

Fp = get_DFT_at_pilots('NoInterf',p,nT);
F = get_DFT_matrix(p);
FBig = kron(eye(nT),F(:,1:numPaths));


nB = 1;
blockLen = p.M*p.K+p.Ncp;
sGFDM = zeros(nB * blockLen, 1, nT);
Cstamp2 = zeros(length(snr),6);
disp_counter = 0;

A_precode = get_transmitter_matrix_precode(p);
A_old = get_transmitter_matrix(p);

pilot_subcar = p.Kset(1)+1:p.Delta_k:p.Kset(end)+1;
data_subcar = p.Kset+1;
data_subcar(pilot_subcar) = [];
M = p.M; K= p.K;
A = zeros(size(A_old));
for m = 1:M
    A(:,data_subcar+K*(m-1)) = A_old(:,data_subcar+K*(m-1)); 
    A(:,pilot_subcar+K*(m-1))= A_precode(:,pilot_subcar+K*(m-1));
end


for si = 1:1%length(snr)
    for i = 1:its
        % Initializations
        s = zeros(p.Kon*p.M-p.Mp*(p.Kon/Delta_k),nT); 
        dd = zeros(p.Kon*p.M-p.Mp*(p.Kon/Delta_k),nT);
        Dd = zeros(p.K,p.M,nT); Dp = zeros(p.K,p.M,nT);
        D = zeros(p.K,p.M,nT);
        x_precode = zeros(N,nT); 
        x_gfdm = zeros(N,nT);
        x_ofdm = zeros(N,nT);
        xp_iT = zeros(N,nT); Xp_iT = zeros(p.K/Delta_k,nT);
        X = zeros(N,nT*N);
        
        for blk = 1:nB
            for iT = 1:nT
                
                % create data symbols
                s(:,iT) = get_random_symbols_new(p); 
                % QAM modulation
                dd(:,iT)= do_qammodulate(s(:,iT), p.mu); 
             
                Dd(:,:,iT) = do_map_p(p, dd(:,iT)); % map them to the D matrix
                                
                % Initialize pilots matrix
                [~, Dp(:,:,iT)] = do_pilot_symbols_time(p, Dd(:,:,iT),Delta_k);
                
                % Space-Frequency pilots insertion
                Dp(:,:,iT) = circshift(Dp(:,:,iT),iT-1,2); 
                
                D(:,:,iT) = Dd(:,:,iT) + Dp(:,:,iT);
                
                xp_iT(:,iT) = do_modulate_precode(p, Dp(:,:,iT),'F', 1:nT);  
                
                Xp_iT(:,iT) = Fp(:,:,iT)*(xp_iT(:,iT)); 
                
                x_precode(:,iT) = do_modulate_precode(p, D(:,:,iT),'F', 1:nT);
                
                x_gfdm(:,iT) = do_modulate(p, D(:,:,iT));
                
                D_temp = D(:,:,iT).';
                D_ofdm = D_temp(:);
                x_ofdm(:,iT) = do_modulate(p_ofdm,D_ofdm);
                
                
            end
            
        end
        
        signal_precode = x_precode(:);
        PAPR_precode(i) = max(abs(signal_precode).^2)/mean(abs(signal_precode).^2);
        
        signal_gfdm = x_gfdm(:);
        PAPR_gfdm(i) = max(abs(signal_gfdm).^2)/mean(abs(signal_gfdm).^2);
        
        signal_ofdm = x_ofdm(:);
        PAPR_ofdm(i) = max(abs(signal_ofdm).^2)/mean(abs(signal_ofdm).^2);
        
        signal_gaussian = randn(N,1);
        PAPR_gaussian(i) = max(abs(signal_gaussian).^2)/mean(abs(signal_gaussian).^2);
        
    
        warning('off','last')
        
    end
    
    % Display some usefull info
    if(mod(disp_counter,4)==0)
        disp_para = sprintf('%dx%d MIMO, K = %d, M = %d, L = %d, alpha = %0.1f',nT,nR,p.K,p.M,numPaths,p.a);  
        disp(disp_para)
        disp_counter = disp_counter + 1; 
    else
        disp_counter = disp_counter + 1; 
    end    
    Cstamp2(si,:) = clock;
    texttodisp = sprintf('SNR = %0.2f with %d iterations ended at %02d:%02d',snr(si),its,Cstamp2(si,4),Cstamp2(si,5));
    disp(texttodisp)

end

plot_dist_l = 1;
plot_dist_s = 1; 
plot_set = [1:plot_dist_l:its*.95, its*.95:plot_dist_s:its]; size(plot_set)
[CDF_precode, x_cdf_precode] = ecdf(PAPR_precode);
figure; semilogy(10*log10(x_cdf_precode(plot_set)),1-CDF_precode(plot_set))

[CDF_gfdm, x_cdf_gfdm] = ecdf(PAPR_gfdm);
if length(x_cdf_gfdm) < its
 CDF_gfdm = [CDF_gfdm; zeros(its - length(x_cdf_gfdm),1)];
 x_cdf_gfdm = [x_cdf_gfdm; zeros(its - length(x_cdf_gfdm),1)];
end
hold on; semilogy(10*log10(x_cdf_gfdm(plot_set)),1-CDF_gfdm(plot_set),'r')

[CDF_gaussian, x_cdf_gaussian] = ecdf(PAPR_gaussian);
hold on; semilogy(10*log10(x_cdf_gaussian(plot_set)),1-CDF_gaussian(plot_set),'k')

[CDF_ofdm, x_cdf_ofdm] = ecdf(PAPR_ofdm);
if length(x_cdf_ofdm) < its
 CDF_ofdm = [CDF_ofdm; zeros(its - length(x_cdf_ofdm),1)];
 x_cdf_ofdm = [x_cdf_ofdm; zeros(its - length(x_cdf_ofdm),1)];
end
hold on; semilogy(10*log10(x_cdf_ofdm(plot_set)),1-CDF_ofdm(plot_set),'m')


legend('Precoded','GFDM','Gaussian','OFDM','Location','southwest')
xlabel('PAPR [dB]')
ylabel('CCDF')
grid





