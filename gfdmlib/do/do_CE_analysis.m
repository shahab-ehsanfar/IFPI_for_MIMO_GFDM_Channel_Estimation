function [ error1, error2, error1z, error2z] = do_CE_analysis(p, Dd, Dp, H, Hhat1, Hhat2, numPaths, snr, its, y)
% Performs analysis of the Channel Estimation in GFDM
%
% Dd: Data Matrix with zero on pilot positions
% Dp: Pilot Matrix with zero on data positions
% H: Channel Frequency response
% Hhat1: Estimated Channel due to 1 pilot
% Hhat2: Estimated Channel due to 2 pilots (Harmonic Mean)
%
% Author: Shahab Ehsanfar




xp1 =  do_modulate(p, [Dp(:,1:p.M-1) Dp(:,p.M-1)]);
xp2 =  do_modulate(p, [Dp(:,2) Dp(:,2:p.M)]);

Xp1 = fft(xp1);
Xp2 = fft(xp2);
% Xp1 = fft(do_modulate(p, Dp));
% Xp2 = Xp1;

xd = do_modulate(p, Dd);
Xd = fft(xd);



%% Harmonic Mean

% Interference Term
InterF = 0.5*H.*(Xd+Xp2)./Xp1 + 0.5*H.*(Xd+Xp1)./Xp2;
interf = ifft(InterF);
interf(numPaths+1:end)=0;
InterFh = fft(interf);

% (Enhanced) Noise term
Noise = Hhat2 - InterFh -H;

if 0 %((snr==5 || snr==30) && its == 1)
    % Comparisons
    figure; semilogy(abs(H),'r')
    hold on; semilogy(abs(H+InterFh+Noise),'k')
    semilogy(abs(InterFh))
    semilogy(abs(Noise),'m')
    semilogy(abs(InterFh+Noise),'g')
    lh = legend('$H$','$\hat{H}$','$\hat{I}_{\rm HM}$', '$\hat{Z}_{\rm HM}$','$\hat{I}_{\rm HM} + \hat{Z}_{\rm HM}$');
    set(lh, 'Interpreter', 'latex')
    xlabel('Frequency (Hz)'); ylabel('$|\cdot|$','interpret','latex');
    string = sprintf('Harmonic Mean, SNR = %d',snr);
    title(string)
    grid
end

%% ONE pilot subsymbol per subcarrier

% Interference term
InterF1 = H.*(Xd./Xp1);
interf1 = ifft(InterF1);
interf1(numPaths+1:end)=0;
InterFh1 = fft(interf1);

% Noise term
Noise1 = Hhat1 - InterFh1 -H;

if 0 % ((snr==5 || snr==30)&& its == 1)
    % Comparisons
    figure;
    semilogy(abs(H),'r')
    hold on; semilogy(abs(H+InterFh1+Noise1),'k')
    semilogy(abs(InterFh1))
    semilogy(abs(Noise1),'m')
    semilogy(abs(InterFh1+Noise1),'g')
    lh = legend('$H$','$\hat{H}$','$\hat{I}_{\rm OP}$', '$\hat{Z}_{\rm OP}$','$\hat{I}_{\rm OP} + \hat{Z}_{\rm OP}$');
    set(lh, 'Interpreter', 'latex')
    xlabel('Frequency (Hz)'); ylabel('$|\cdot|$','interpret','latex');
    string = sprintf('One Pilot, SNR = %d',snr);
    title(string)
    grid
end

%% Compare Harmonic Mean vs One Pilot
if 0 %((snr==5 || snr==30)&& its == 1)
figure; semilogy(abs(InterFh1+Noise1).^2)
hold on; semilogy(abs(InterFh+Noise).^2,'r')
legend('One Pilot','Harmonic Mean')
xlabel('Frequency (Hz)'); ylabel('$|\hat{I} + \hat{Z}|^2$','interpret','latex');
string = sprintf('SNR = %d',snr);
title(string)
grid
end

%% Alternative MSE computation for a single run
error1 = sum(abs(InterFh1).^2)./(p.M*p.K);
error2 = sum(abs(InterFh).^2)./(p.M*p.K);

error1z = sum(abs(Noise1).^2)./(p.M*p.K);
error2z = sum(abs(Noise).^2)./(p.M*p.K);


%% Plot Impulse response of Noise Enhancements

% figure; plot(abs(ifft(1./Xp1)))
% figure; plot(abs(ifft(0.5*(Xp1+Xp2)./(Xp1.*Xp2))))


end

