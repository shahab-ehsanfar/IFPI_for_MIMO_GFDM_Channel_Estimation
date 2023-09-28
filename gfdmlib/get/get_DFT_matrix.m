function F = get_DFT_matrix(p)

N= p.M*p.K;

F = zeros(N,N);
for k = 0:N-1
    for i = 0:N-1
        F(i+1,k+1) = exp(-2j*pi*k*i/N);
    end
end

F = (1/sqrt(N))*F;
