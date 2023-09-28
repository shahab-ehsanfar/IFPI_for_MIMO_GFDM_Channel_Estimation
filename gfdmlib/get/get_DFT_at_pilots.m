function Fp = get_DFT_at_pilots(mode,p,nT)

if nargin < 3
    nT = 1;
end
F = get_DFT_matrix(p);
Delta_k = p.Delta_k;
switch mode
    case 'NoInterf'
        
        kset = get_kset(p);       
        
        Fp_i = zeros(length(kset)/Delta_k,p.M*p.K);
        Fp = zeros(length(kset)/Delta_k,p.M*p.K,nT);
        
        for iT = 1:nT
            for k = 0:length(kset)/Delta_k-1
                Fp_i(k+1,:) = F(1+kset(1)*p.M+(iT-1)+k*Delta_k*p.M,:);
            end
            Fp(:,:,iT) = Fp_i(:,:);
        end
        
        
    case 'exact'
        for k = 0:p.K/Delta_k
            if k == 0
                for m = 1:p.M-floor(p.M/2)
                    Fp(m,:) = F(m,:);
                end
            elseif k == p.K/Delta_k
                for m = 1:p.M-ceil(p.M/2)
                    i_out = m+k*Delta_k*p.M-floor(p.M/2);
                    i_in = m+k*p.M-floor(p.M/2);
                    Fp(i_in,:) = F(i_out,:);
                end
            else
                for m = 1:p.M
                    i_out = m+k*Delta_k*p.M-floor(p.M/2);
                    i_in = m+k*p.M-floor(p.M/2);
                    Fp(i_in,:) = F(i_out,:);
                end
            end
        end
        
    case 'wider'       
        
        Fp = zeros(p.M*p.K/Delta_k,p.M*p.K);
        for k = 0:p.K/Delta_k
            if k == 0
                for m = 1:p.M-floor(p.M/2)+1
                    Fp(m,:) = F(m,:);
                end
            elseif k == p.K/Delta_k
                for m = 0:p.M-ceil(p.M/2)
                    i_out = m+k*Delta_k*p.M-floor(p.M/2);
                    i_in = m+1+k*(p.M+2)-floor((p.M+2)/2);
                    Fp(i_in,:) = F(i_out,:);
                end
            else
                for m = 0:p.M+1
                    i_out = m+k*Delta_k*p.M-floor(p.M/2);
                    i_in = m+1+k*(p.M+2)-floor((p.M+2)/2);
                    %[i_out i_in]
                    Fp(i_in,:) = F(i_out,:);
                end
            end
        end
        
end