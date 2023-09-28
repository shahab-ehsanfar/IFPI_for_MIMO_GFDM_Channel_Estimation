% Copyright (c) 2014 TU Dresden
% All rights reserved.
% See accompanying license.txt for details.
%
clear variables

snr = -5:5:35;
p=get_defaultGFDM('OFDM');
p.pulse='rc_fd';
p.a=0.3;%roll-off factor for the RRC and RC filters.
p.K=96;%Number of useful subcarriers.
p.L=2;%Number of overlapping subcarrier for frequency domain GFDM implementation..
p.Kset = 0:p.K(end)-1;
p.Kon = length(p.Kset);
%p.Kon=floor(p.K/2)+1;%Number of actives subcarriers.
p.Mp = 2; % Number of pilot subsymbols
p.Md = 5; % Number of data subsymbols 
p.M= p.Md + p.Mp; % Number of subsymbols.
p.mu=4; % 16-QAM.
p.Ncp = p.K; % Number of CP samples
p.b = 0; %p.K/4; % Number of window raise/fall samples
p.window = 'rc';
%p.window = 'rc';
p.mode = 'precode';

its = 30; % Simulation iterations for each SNR

% Currently pilot insertion is adapted for only up to 2x2 MIMO channel
nT = 2; % Number of Tx antennas (up to 2)
nR = 2; % Number of Rx antennas (up to 2)

p.Delta_k = 3; % pilot subcarrier spacing
numPaths = 24; % number of non-zero channel taps

N = p.M*p.K;
Delta_k = p.Delta_k;

errorsZF_gen = zeros(length(snr), its);
errorsZF_ls = zeros(length(snr), its);
errorsZF_lmmse = zeros(length(snr), its);


P = get_power_delay_profile(numPaths,nT,nR); 

Fp = get_DFT_at_pilots('NoInterf',p,nT);
F = get_DFT_matrix(p);
FBig = kron(eye(nT),F(:,1:numPaths));

% Use the same channel auto correlation for all antennas
% as PDP is assumed the same for all antennas
P_freq = sqrt(N)*Fp(:,:,1)*diag([sqrt(P(:,1,1)); zeros(N-numPaths,1)]);
R_HH = P_freq*P_freq';

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
    A(:,data_subcar+K*(m-1)) = A_old(:,data_subcar+K*(m-1)); %*(sqrt(M)/sqrt(M-1));
    A(:,pilot_subcar+K*(m-1))= A_precode(:,pilot_subcar+K*(m-1));
end


for si = 1:length(snr)
    for i = 1:its
        % Initializations
        s = zeros(p.Kon*p.M-p.Mp*(p.Kon/Delta_k),nT);  
        dd = zeros(p.Kon*p.M-p.Mp*(p.Kon/Delta_k),nT);
        Dd = zeros(p.K,p.M,nT); Dp = zeros(p.K,p.M,nT);
        D = zeros(p.K,p.M,nT);
        x_iT_noCP = zeros(N,nT); X_iT = zeros(N,nT);
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
                Dp(:,:,iT) = circshift(Dp(:,:,iT),iT-1,2); %*2.3535;
                
                D(:,:,iT) = Dd(:,:,iT) + Dp(:,:,iT);
                
                xp_iT(:,iT) = do_modulate_precode(p, Dp(:,:,iT),'F', 1:nT);  
                
                Xp_iT(:,iT) = Fp(:,:,iT)*(xp_iT(:,iT)); 
                
                x_iT_noCP(:,iT) = do_modulate_precode(p, D(:,:,iT),'F', 1:nT);
                
                X_iT(:,iT) = fft_u(x_iT_noCP(:,iT));
                
                
                % X in wide matrix form                
                for n = 1:N
                    X(n,N*(iT-1)+n) = X_iT(n,iT);
                end

                
            end
            
        end
        
        % Random channel
        h = get_mimo_channel(P,0,numPaths,0,nT,nR); 
       
        % AWGN vector in frequency
        W = fft_u(awgn(zeros(N,nR)*(0+0i),snr(si)));
        
        % Receive signal in frequency
        Y = sqrt(N)*X*FBig*h + W;
        
        [Heq_ls, H_ls] = do_Chan_Estimation_mimo_precode('LS', p, Y, numPaths, snr(si), N, nT, nR, Fp, Xp_iT, h);
        
        [Heq_lmmse, H_lmmse] = do_Chan_Estimation_mimo_precode('LMMSE', p, Y, numPaths, snr(si), N, nT, nR, Fp, Xp_iT, h, R_HH);
        
        Heq = zeros(N*nT,N*nR);
        H_i = zeros(N,nT,nR);
        for iR = 1:nR
            for iT = 1:nT
                H_i(:,iT,iR) = fft(h((iT-1)*numPaths+1:(iT-1)*numPaths+numPaths,iR),N);
                for n = 1:N
                    Heq((iR-1)*N+n,(iT-1)*N+n) = H_i(n,iT,iR);
                end
            end
        end
        
        mse_ls(i,si) = mean(abs(H_i(:)-H_ls(:)).^2);
        mse_lmmse(i,si) = mean(abs(H_i(:)-H_lmmse(:)).^2);
        
     
        % Channel equalization for 2x2 MIMO
        [Xhat_genie,~] = do_chan_equalization(N,Heq,Y,snr(si),nT);
        [Xhat_ls,~] = do_chan_equalization(N,Heq_ls,Y,snr(si),nT);        
        [Xhat_lmmse,~] = do_chan_equalization(N,Heq_lmmse,Y,snr(si),nT);        
        
        Dhat_genie = zeros(p.K,p.M,nT);
        Dhat_ls = Dhat_genie; Dhat_lmmse = Dhat_genie;
        dhat_genie = zeros(p.Kon*p.M-p.Mp*(p.Kon/Delta_k),nT);
        dhat_ls = dhat_genie; dhat_lmmse = dhat_genie;
        sh_genie = zeros(p.Kon*p.M-p.Mp*(p.Kon/Delta_k),nT); 
        sh_ls = sh_genie; sh_lmmse = sh_genie;
        for iT =1:nT
            % Zero Forcing Receiver
            Dhat_genie(:,:,iT) = do_demodulate_combine(p, ifft_u(Xhat_genie(N*(iT-1)+1:iT*N)),'ZF',1:nT);                     
            Dhat_ls(:,:,iT) = do_demodulate_combine(p, ifft_u(Xhat_ls(N*(iT-1)+1:iT*N)),'ZF',1:nT);         
            Dhat_lmmse(:,:,iT) = do_demodulate_combine(p, ifft_u(Xhat_lmmse(N*(iT-1)+1:iT*N)),'ZF',1:nT);         
                        
            % Unmap
            dhat_genie(:,iT) = do_unmap_p(p, Dhat_genie(:,:,iT));
            dhat_ls(:,iT) = do_unmap_p(p, Dhat_ls(:,:,iT));
            dhat_lmmse(:,iT) = do_unmap_p(p, Dhat_lmmse(:,:,iT));  
            
            % QAM demodulate
            sh_genie(:,iT) = do_qamdemodulate(dhat_genie(:,iT),p.mu);
            sh_ls(:,iT) = do_qamdemodulate(dhat_ls(:,iT),p.mu);
            sh_lmmse(:,iT) = do_qamdemodulate(dhat_lmmse(:,iT),p.mu);
        end     
        
        errorsZF_gen(si,i) = sum(sh_genie(:) ~= s(:));       
        errorsZF_ls(si,i) = sum(sh_ls(:) ~= s(:));            
        errorsZF_lmmse(si,i) = sum(sh_lmmse(:) ~= s(:));            
        
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

Es = 10*log10(1+(p.Mp-1)/(p.Md+(Delta_k-1)*p.M)  );

serZF_gen = sum(errorsZF_gen.');
serZF_ls = sum(errorsZF_ls.');
serZF_lmmse = sum(errorsZF_lmmse.');


% Plot the results

figure;
semilogy(snr+Es,mean(mse_ls),'b-s','LineWidth',1.5)
hold on;
semilogy(snr+Es,mean(mse_lmmse),'r-','LineWidth',1.5)
grid
xlabel('SNR (dB)')
ylabel('MSE')
legend({'LS', 'LMMSE'});

figure
semilogy(snr, serZF_gen./(its*p.M*p.K*nT),'b');
hold on;
semilogy(snr+Es, serZF_ls./(its*p.M*p.K*nT),'r');
semilogy(snr+Es, serZF_lmmse./(its*p.M*p.K*nT),'m');
grid;
xlabel('SNR (dB)')
ylabel('SER')
legend({'Genie Aided', 'LS', 'LMMSE'});





