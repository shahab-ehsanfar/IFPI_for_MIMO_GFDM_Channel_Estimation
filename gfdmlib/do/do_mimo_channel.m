% Copyright (c) 2014 TU Dresden
% All rights reserved.
% See accompanying license.txt for details.
%


function xch = do_mimo_channel(mode, x_mimo, h, numPaths, nT, nR, signal_length, noblk)
% Apply channel convolution and noise to a signal X
%
% x_mimo: transmit signal
% h: matrix of channel imp. resp.
% nT: number of transmit antennas
% noblk: number of transmitted gfdm data blocks
% signal_length: including cp length
%

switch mode
    case 'Filter'
        xch = zeros(noblk*signal_length,nR);
        xch_iRiT = zeros(noblk*signal_length,nT);        
        for iR = 1:nR
            for iT = 1:nT
                start_pnt = (iT-1)*numPaths+1;
                end_pnt = start_pnt + numPaths-1;
                                
                h_ii = h(start_pnt:end_pnt,iR);
                
                % Filter the signal with channel imp. resp.
                xch_iRiT(:,iT) = filter(h_ii,1,x_mimo(:,iT));
            end
            xch(:,iR) = sum(xch_iRiT,2);
        end
        
    case 'Filter time variant'
        xch = zeros(noblk*signal_length,nR);
        xch_iRiT = zeros(noblk*signal_length,nT);
        for iR = 1:nR
            for iT = 1:nT                                
                % Filter the signal with channel imp. resp.
                xch_iRiT(:,iT) = filter(h(iT,iR),x_mimo(:,iT));
            end
            xch(:,iR) = sum(xch_iRiT,2);
        end
        
    case 'Toeplitz' % Slower because of matrix multiplication
        
        h_matrix = [];
        for iR = 1:nR
            h_matrix_iR = [];
            for iT = 1:nT
                start_pnt = (iT-1)*numPaths+1;
                end_pnt = start_pnt + numPaths-1;
            
                h_long = ifft(fft(h(start_pnt:end_pnt,iR),noblk*signal_length));
                
                % Create lower triangular channel matrix in time
                h_matrix_iRiT = toeplitz(h_long,[h_long(1); zeros(length(h_long)-1,1)]);                 
                h_matrix_iRiT(abs(h_matrix_iRiT) < 10e-16) = 0;
                
                % construct the mimo channel matrix in time
                h_matrix_iR = [h_matrix_iR h_matrix_iRiT];
            end
            h_matrix = [h_matrix; h_matrix_iR];
        end        
        xch = h_matrix*x_mimo(:);
    
end


















        
        
        
       
        