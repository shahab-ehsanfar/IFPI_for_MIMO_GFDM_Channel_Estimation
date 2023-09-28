function [ C_h,E_h ] = calc_softdec_equichan(p,Ainv_F_mimo,F_mimo_A,QHeq,Q,nT)


N = p.K*p.M;        
if (nT == 2 && strcmp(p.mode,'precode'))  % Efficient 2x2 calculations 
        AFQ = zeros(size(Q)); 
        AFQH = zeros(size(Q));        
        for n = 1:N
            AFQ(:,[n n+N]) = Ainv_F_mimo(:,[n n+N])*Q([n n+N],[n n+N]);
            AFQH(:,[n n+N]) = Ainv_F_mimo(:,[n n+N])*QHeq([n n+N],[n n+N]);
        end
        Sigma_genie = sum(AFQ.*conj(AFQ),2);
      
        F_mimo_A_spy = F_mimo_A;
        F_mimo_A_spy(abs(F_mimo_A_spy) < 1e-10) = 0;
        AFQH_spy = AFQH;
        AFQH_spy(abs(AFQH_spy) < 1e-10) = 0;
        E = full( sparse(AFQH_spy)*sparse(F_mimo_A_spy) );
        Ev = diag(E);
        
elseif strcmp(p.mode,'ofdm')
    diagQHeq = Q;
    Es =1 ; % mean(dd(:,1) .* conj(dd(:,1)) + dd(:,2) .* conj(dd(:,2)) )/2; %  mean(dd(:,1) .* conj(dd(:,1)));
    Sigma_genie = diagQHeq - Es ; %1; 
    
    %E = QHeq;
    
    Ev = ones(nT*N,1);%diag(E).*diagQHeq;
   
    
else % Exact matrix calculation
    Sigma_genie = diag(Ainv_F_mimo*(Q*Q')*Ainv_F_mimo'); % SLOWEST LINE
    E = (Ainv_F_mimo*QHeq*F_mimo_A);
    Ev = diag(E);
    

end
        
        siSigma_genie = 1./sqrt(Sigma_genie);
        
        
        d_set = 1:p.M*p.Kon*nT;
        if p.Mp == 1
            d_set(1:p.Delta_k:p.Kon) = [];
        elseif (p.Mp == 2 && nT == 2 && strcmp(p.mode,'precode'))
            d_set([1:p.Delta_k:p.Kon p.Kon+1:p.Delta_k:p.Kon+p.Kon...
                    N+(1:p.Delta_k:p.Kon) N+(p.Kon+1:p.Delta_k:p.Kon+p.Kon) ]) = [];            
        elseif (p.Mp == 2 && nT == 1 && strcmp(p.mode,'precode'))
            d_set([1:p.Delta_k:p.Kon p.Kon+1:p.Delta_k:p.Kon+p.Kon]) = []; 
        elseif p.Mp == 2 && strcmp(p.mode,'ofdm')
            
           % d_set([1:p.Delta_k*p.M:N*nT 2:p.Delta_k*p.M:N*nT]) = [];
            %    d_set = 1:96*7; d_set([1:14:96*7 14:14:96*7]) = [];     
        else
            error('PLEASE CHECK DATA SET d_set')
        end
        
        C_h = siSigma_genie(d_set).';
        E_h = Ev(d_set).';
 
end

