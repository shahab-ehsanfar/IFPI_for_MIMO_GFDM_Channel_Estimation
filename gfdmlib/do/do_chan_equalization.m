function [ Xhat, Q, QHeq ] = do_chan_equalization(N,Heq,Y,snr,nT,p)
% Works only for 2x2 MIMO

Xhat = zeros(size(Heq,2), 1);
Y_long = Y(:);
 VarN = 1/(10^(snr/10));
Q = zeros(size(Heq));
QHeq = zeros(size(Heq));

if nargin < 6
    if nT == 2
        % Efficient 2x2 equalization
        for n =1:N
            % Equalization
            Hpart = Heq([n, n+N], [n, n+N]);
            Ypart = Y_long([n, n+N]);
            invPart = (Hpart'*Hpart+VarN*eye(2));
            filter = invPart \ Hpart';
            Xhat([n, n+N]) = (filter*Ypart); %./diag(filter);
            
            % Extra information for soft decoding:
            Q([n, n+N], [n, n+N]) = filter; %(filter*filter'); %./diag(diag(filter).^2);
            QHeq([n, n+N], [n, n+N]) = filter*Hpart;
        end
    else
        % General form
        Q = (Heq'*Heq+VarN*eye(nT*N)) \ Heq'; %diag(1./diag(Heq));%
        
        QHeq = Q*Heq;
        
        Xhat = Q*Y(:); %./diag(filter2);
    end
else
        
    if nT == 2
        diagQHeq = zeros(size(Xhat));
        for n =1:N
            % Equalization
            Hpart = Heq([n, n+N], [n, n+N]);
            Ypart = Y_long([n, n+N]);
            invPart = (Hpart*Hpart'+VarN*eye(2));
            filter = Hpart'/invPart;
            diagQHeq([n, n+N]) = 1./diag(filter*Hpart);
            
            Xhat([n, n+N]) = (filter*Ypart).*diagQHeq([n, n+N]); %./diag(filter);
            
            % Extra information for soft decoding:
            %Q([n, n+N], [n, n+N]) = filter; %(filter*filter'); %./diag(diag(filter).^2);
            QHeq([n, n+N], [n, n+N]) = filter*Hpart;
        end
    else
        Q = Heq'/(Heq*Heq'+ VarN*eye(nT*N));
        
        diagQHeq = 1./diag( Q*Heq );
    
        Xhat = (Q*Y(:)).*diagQHeq;
        
        QHeq = (Q*Heq);
    end
    
    
    Q = diagQHeq;
       
    
end
 
end

