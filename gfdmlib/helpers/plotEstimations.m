function [ PLOT ] = plotEstimations( H, Hp_mf, Hp_zf, numPaths, snr, numplot)

switch numplot
    case 3
        PLOT = figure;
        plot(10*log10(abs(H)),'r','LineWidth',2)
        hold on
        plot(10*log10(abs(Hp_mf)),'--','LineWidth',2)
        plot(10*log10(abs(Hp_zf)),'g-.','LineWidth',2)
        legend('Channel','1 pilot Estimation','2 pilots Estimation')
        string = sprintf('SNR = %d dB, no. of paths: %d',snr,numPaths);
        title(string)
        xlabel('Frequency (Hz)')
        ylabel('|H_p| (dB)')
    case 2
        PLOT = figure;
        semilogy(abs(H),'r','LineWidth',1)
        hold on
        semilogy(abs(Hp_mf),'LineWidth',1)
        legend('Channel','Estimation')
        string = sprintf('SNR = %d dB, no. of paths: %d',snr,numPaths);
        title(string)
        xlabel('Frequency (Hz)')
        ylabel('|H_p| (dB)')
end

end

