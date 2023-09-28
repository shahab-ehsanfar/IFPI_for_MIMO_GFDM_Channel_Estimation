function [Hp]= do_channelEstimation(p, Dx, Dy, method, dim)
% Estimates the channel based on pilot symbols
%
% Dx is the matrix of data and pilot symbols before modulation
% Dy is the matrix of received data and pilot symbols
%
% Author: Shahab Ehsanfar

switch dim
    case 'freq'
        % Extract the pilot subcarriers
        ksetp = (p.pf:p.pf:p.K);
        
        Xp = fft(Dx(ksetp,:),[],2);
        Yp = fft(Dy(ksetp,:),[],2);
        
        % Channel Estimation
        switch method
            case 'LS' % Least Square
                Hp_i = Yp ./ Xp;
                
                %     case 'LMMSE' % Linear Minimum Mean Square Error
                %         Hp_ls = Yp ./ Xp;
                %
                %         % channel impulse response based on pilot subcarriers
                %         hp_ls = ifft(Hp_ls,[],2);
                %
                %         % Average power of each tap for p.M subsymbols
                %         P_ls = (1/p.M)*sum( abs(hp_ls.').^2 );
                %
                %         % Find the most significant taps
                %         [~ ,index] = sort(P_ls,2);
                %         MST = index( end-round(p.Kpilot/4) : end)'; % The 1/4 pilot subcarriers with the largest powers
                %
                %         % Select the 1/4 most significant taps
                %         P_mst = zeros(p.Kpilot,1);
                %         P_mst(MST) = P_ls(MST);
                %
                % %         % The channel autocorrelation matrix is obtained by circular shift
                % %         % of A_tilde
                % %         A_tilde = p.Kpilot * ifft(P_mst);
                % %         for k = 1:p.Kpilot
                % %             R_hphp(k,:) = circshift(A_tilde.',k-1,2);
                % %         end
                %
                %         % Approximate SNR
                %         SNR_tilde = sum(P_mst)/(sum(P_ls)-sum(P_mst));
                %
                %         switch p.mu
                %             case 2
                %                 beta = 1;
                %             case 4
                %                 beta = 17/9;
                %             case 6
                %                 beta = 2.6854;
                %         end
                %
                %         % The first row B_tilde of the LMMSE matrix
                %         row = zeros(1,p.Kpilot);
                %         for k = 1:p.Kpilot
                %            row(k) = P_mst(k)/( P_mst(k) + (beta/(p.Kpilot*SNR_tilde)) );
                %         end
                %         B_tilde = ifft(row,[],2);
                %
                %         % LMSSE is obtained by circular shift of B_tilde
                %         lmmse = zeros(p.Kpilot,p.Kpilot);
                %         for k = 1:p.Kpilot
                %             lmmse(k,:) = circshift(B_tilde,k-1,2);
                %         end
                %
                % %         lmmse = R_hphp* inv(R_hphp + (beta/SNR_tilde)*eye(p.Kpilot));
                %
                %         Hp_i = lmmse*Hp_ls;
        end
        Hp_mean = (1/p.M)*sum( (Hp_i.') );
        
        % copy the last element to the begenning in order to make it look periodic
        Hp_mean= [Hp_mean(end) Hp_mean];
        
        % Perform the channel interpolation
        int_o = (1:p.K*p.M);
        int_i = [1,p.pf*p.M : p.M*p.pf : p.K*p.M];
        Hp_interp = pchip(int_i,Hp_mean,int_o).'; % Piecewise Cubic Hermite Interpolating Polynomial
        %Hp_interp = interp1(int_i,Hp_mean,int_o).';
        Hp = circshift(Hp_interp, -p.M+1);
    
    
    case 'time'
        
        if p.Mp == 2        
            msetp = [1 p.M];
        else
            msetp = 1;
        end
        
        % Perform FFT and then Extract the pilot subcarriers  
        X = fft(Dx,[],2);
        Y = fft(Dy,[],2);        
        Xp = X(:,msetp);
        Yp = Y(:,msetp);
              
        % Channel Estimation
        switch method
            case 'LS' % Least Square
                Hp_i = Yp ./ Xp;
        end
        
        if p.Mp == 2
            Hp_mean = (1/p.Mp)*sum( (Hp_i.') );
        else
            Hp_mean = Hp_i.';
        end
        
        % copy the last element to the begenning in order to make it look periodic
        Hp_mean= [Hp_mean(end) Hp_mean];
        
        % Perform the channel interpolation
        int_o = (1:p.K*p.M);
        int_i = [1,p.M:p.M: p.K*p.M];
        % type 'pchip' for Piecewise Cubic Hermite Interpolating Polynomial 
        % type 'linear' for linear interpolation
        Hp_interp = interp1(int_i,Hp_mean,int_o,'pchip').';  
        Hp = circshift(Hp_interp, -p.M+1);
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        

end

       
end

