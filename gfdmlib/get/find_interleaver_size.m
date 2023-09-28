function [PCCC, interleaver, bicm] = find_interleaver_size(p, Nd_blk, nT)
% Searches the possible interleavers for the given number of transmit
% symbols
% Nd_blk: Number of QAM symbols per block for information transfer
% nT: number of Tx antennas
% p: GFDM structure
% This is a beta version (still may need to be debugged)
% Shahab 

Pathadd = pwd;
if  strcmp(Pathadd(1),'/')
    PCCC = load('./Lib/cml/PCCCk88.mat');
else
    %PCCC = load('C:\Documents\MATLAB\GFDM\Lib\cml\PCCCk88.mat');
    PCCC = load('C:\Users\Shahab Ehsan Far\Documents\MATLAB\GFDM\Lib\cml\PCCCk88.mat');
end

code_param = PCCC.code_param;
sim_param = PCCC.sim_param;

% create S-matrix for 2^mu-QAM
x = (0:(2^p.mu-1)).';
S_matrix0 = qammod(x,2^p.mu).';
S_matrix = S_matrix0/sqrt( mean( abs(S_matrix0).^2 ) ); % Normalize to unit power
code_param.S_matrix = S_matrix;


Nd = Nd_blk*nT*p.mu;
N_tail = -2;
N_tail_old = N_tail;
i_try = 0;
code_length =0;

if p.coderate == 3/4
    if nT == 2
      % code_param.pun_pattern = [1,0,1,0,1,0;1,0,1,0,1,0;0,0,0,0,0,0;1,0,1,0,0,0]; 
        code_param.pun_pattern = [1,0,1,0,1,0;1,0,0,1,0,0;0,0,1,0,0,1;0,0,0,1,0,0];
    else
        code_param.pun_pattern = [1,0,1,0,1,0;1,0,0,1,0,0;0,0,1,0,0,1;0,0,0,1,0,0];%[1,0,1,0,1,0;1,0,1,0,1,0;0,0,0,0,0,0;1,0,1,0,0,0]; 
    end
elseif p.coderate == 2/3
   code_param.pun_pattern = [1,0,1,0;0,0,0,1;0,1,0,0;0,0,1,1];   
end

code_param.symbols_per_frame = Nd_blk*nT;
code_param.code_bits_per_frame = Nd;


if (p.coderate == 1/3 || p.coderate == 3/4 || p.coderate == 2/3)
    while(code_length ~= Nd)
        Nb_test = round(Nd*p.coderate) + N_tail;
        data = round(rand([1 Nb_test]));
        code_interleaver = randperm(Nb_test)-1;
        codeword = TurboEncode( data, code_interleaver, code_param.pun_pattern, code_param.tail_pattern, sim_param.g1, sim_param.nsc_flag1, sim_param.g2, sim_param.nsc_flag2 );
        code_length = length(codeword);
        if code_length > Nd
            if (code_length - 1) == Nd
                code_param.tail_pattern(end,end) = 0;
            elseif (code_length - 2) == Nd
                code_param.tail_pattern(end-1:end,end) = 0;
            elseif N_tail_old ~= N_tail -1
                N_tail = N_tail - 1;
            end
            
        elseif code_length < Nd
            N_tail_old = N_tail;
            N_tail = N_tail + 1;
            
        end
        if (i_try > 30)
            error('Cannot find Nb after 30 trials')
        end
        i_try = i_try +1;
    end
    
    interleaver = Nb_test ;
    PCCC.code_param = code_param;
    bicm = Nd;
end
