function [ Matrix_out ] = get_blkdiag( Matrix_in )
% Generates block diagonal matrix from a 3D array
%   
        N = size(Matrix_in,3);
        Mat_str = char(zeros(N,15 + numel(num2str(N)) ));
        for n = 1:N
            Mat_str(n,1:15+numel(num2str(n))) = sprintf('Matrix_in(:,:,%d)',n);
            
            if n == 1
                Mat_blks = Mat_str(n,1:15+numel(num2str(n)));
            else
                Mat_blks = [Mat_blks ',' Mat_str(n,1:15+numel(num2str(n)))];
            end
        end
        Matrix_out = eval(['blkdiag(' Mat_blks ');']); 
end

