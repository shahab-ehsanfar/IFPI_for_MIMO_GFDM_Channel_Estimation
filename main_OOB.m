% Copyright (c) 2014 TU Dresden
% All rights reserved.
% See accompanying license.txt for details.
%


% This example shows how a random GFDM signal can be
% created. Additionally, an OFDM signal with an equal length is
% created, and the Spectrum of both is compared.

%% Create parameter sets for GFDM and OFDM
gfdm = get_defaultGFDM('TTI');
gfdm.K = 96;
gfdm.Kset = 20:40;  % Only allocate some subcarriers
gfdm.pulse = 'rc_fd';
gfdm.a = 0.3;
gfdm.M = 15;
gfdm.Mon = gfdm.M;
gfdm.window = 'rc';
gfdm.b = gfdm.K/4;
gfdm.Delta_k = 3;
gfdm.mu = 4;
gfdm.Mp = 2;
gfdm.Md = gfdm.M-gfdm.Mp;
gfdm.Kon = length(gfdm.Kset);
gfdm.mode = 'precode';


ofdm = gfdm;
ofdm.pulse = 'rc_td';  % use RC_TD with rolloff 0 to make a
ofdm.a = 0;            % rectangular filter
ofdm.M = 1;
ofdm.Mon = 1;
ofdm.L = ofdm.K;

its = 20;

nB = 3; % Number of GFDM blocks to generate

%% Generate the signals
% Currently only works without CP
assert(~isfield(gfdm, 'Ncp') || gfdm.Ncp == 0);

% Allocate enough space for the signals
blockLen = gfdm.M*gfdm.K;
sGFDM = zeros(nB * blockLen, its);
sOFDM = zeros(size(sGFDM));
sGFDM_ifpi = zeros(size(sGFDM));
sWGFDM = zeros(size(sGFDM));
sWGFDM_ifpi = zeros(size(sGFDM));

for i = 1:its
    for b = 1:nB
        % Create GFDM signal by modulation of random data
        s = get_random_symbols_new(gfdm);
        dd = do_qammodulate(s, gfdm.mu);
        Dd = do_map_p(gfdm, dd);
        [~, Dp] = do_pilot_symbols_time(gfdm, Dd,gfdm.Delta_k);
        Dp = circshift(Dp,2-1,2);
        
        D = Dd+Dp;
        
        x = do_modulate(gfdm, D);
        w = get_window(gfdm);
        x_w = x.*w;
        sGFDM((b-1)*blockLen+(1:blockLen),i) = x;
        sWGFDM((b-1)*blockLen+(1:blockLen),i) = x_w;
        
        
        x_ifpi = do_modulate_precode(gfdm, D,'F',1:2);
        w = get_window(gfdm);
        x_w_ifpi = x_ifpi.*w;
        sGFDM_ifpi((b-1)*blockLen+(1:blockLen),i) = x_ifpi;
        sWGFDM_ifpi((b-1)*blockLen+(1:blockLen),i) = x_w_ifpi;
        
        
        % Create an OFDM signal
        for m = 1:gfdm.M
            D2 = do_map(ofdm, do_qammodulate(get_random_symbols(ofdm), ofdm.mu));
            x2 = do_modulate(ofdm, D2);
            sOFDM((b-1)*blockLen+(m-1)*gfdm.K+(1:gfdm.K),i) = x2;
        end
    end
end
%% Plot the resulting PSD
f = linspace(-gfdm.K/2, gfdm.K/2, 2*length(sGFDM)+1); f = f(1:end-1)';


psd_ofdm_exact = mean((fftshift(abs(fft(sOFDM, 2*length(sOFDM))))).'/2).';
psd_gfdm_ifpi_exact = mean((fftshift(abs(fft(sGFDM_ifpi, 2*length(sGFDM))))).'/2).';
psd_gfdm_exact = mean((fftshift(abs(fft(sGFDM, 2*length(sGFDM))))).'/2).';
psd_wgfdm_ifpi_exact = mean((fftshift(abs(fft(sWGFDM_ifpi, 2*length(sGFDM))))).'/2).';
psd_wgfdm_exact = mean((fftshift(abs(fft(sWGFDM, 2*length(sGFDM))))).'/2).';

L = 10;
psd_ofdm = abs(ifft(fft(psd_ofdm_exact).*fft(ones(1,L),length(psd_ofdm_exact)).')/sqrt(L));
psd_gfdm_ifpi = abs(ifft(fft(psd_gfdm_ifpi_exact).*fft(ones(1,L),length(psd_ofdm_exact)).')/sqrt(L));
psd_gfdm = abs(ifft(fft(psd_gfdm_exact).*fft(ones(1,L),length(psd_ofdm_exact)).')/sqrt(L));
psd_wgfdm_ifpi = abs(ifft(fft(psd_wgfdm_ifpi_exact).*fft(ones(1,L),length(psd_ofdm_exact)).')/sqrt(L));
psd_wgfdm = abs(ifft(fft(psd_wgfdm_exact).*fft(ones(1,L),length(psd_ofdm_exact)).')/sqrt(L));

c = @cmu.colors;

plot_set = 1:12:length(f);
figure;
plot(f(plot_set), -39+mag2db(psd_ofdm(plot_set)), 'b'); 
hold on;
plot(f(plot_set), -39+mag2db(psd_gfdm_ifpi(plot_set)),'g');%,'color',c('La Salle Green'));
plot(f(plot_set), -39+mag2db(psd_gfdm(plot_set)), 'r');
plot(f(plot_set), -39+mag2db(psd_wgfdm_ifpi(plot_set)),'k');
plot(f(plot_set), -39+mag2db(psd_wgfdm(plot_set)), 'm');
hold off;
ylim([-120, 10]);
xlim([-gfdm.K/2 +gfdm.K/2]);
xlabel('f/F'); ylabel('PSD [dB]');
grid()
legend({'OFDM', 'IFPI GFDM','Basic GFDM','W-IFPI GFDM','W-GFDM'},'Location','southeast');





