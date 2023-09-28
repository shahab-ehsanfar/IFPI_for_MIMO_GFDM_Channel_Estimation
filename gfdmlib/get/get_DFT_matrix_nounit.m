function F = get_DFT_matrix_nounit(p)

N= p.M*p.K;

F = zeros(N,N);
for k = 0:N-1
    for i = 0:N-1
        F(i+1,k+1) = exp(-2j*pi*k*i/N);
    end
end

