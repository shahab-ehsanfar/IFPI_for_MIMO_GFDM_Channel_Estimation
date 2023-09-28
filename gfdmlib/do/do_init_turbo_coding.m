% Copyright (c) 2014 TU Dresden
% All rights reserved.
% See accompanying license.txt for details.
%

function [PCCC1, PCCC2, interleaver, bicm] = do_init_turbo_coding(p)

Pathadd = pwd;

% Load a sample structure for Turbo-coding
if  strcmp(Pathadd(1),'/')
    PCCC1 = load('./cml/PCCCk88.mat');
else
    PCCC1 = load('C:\Documents\MATLAB\GFDM\Lib\cml\PCCCk88.mat');
end
PCCC2 = PCCC1;
if p.K*p.M == 672
    if p.Delta_k == 1 && p.coderate == 1/3
        interleaver.K1 = 764;
        bicm.K1 = 576*p.mu;
        interleaver.K2 = 636;
        bicm.K2 = 480*p.mu;
    elseif p.Delta_k == 2 && p.coderate == 1/3
        interleaver.K1 = 828;
        bicm.K1 = 624*p.mu;
        interleaver.K2 = 764;
        bicm.K2 = 576*p.mu;
    elseif p.Delta_k == 3 && p.coderate == 1/3
        interleaver.K1 = 860;
        bicm.K1 = 640*p.mu;
        interleaver.K2 = 807;
        bicm.K2 = 608*p.mu;
        PCCC2.code_param.tail_pattern(end,end)=0;
    elseif p.Delta_k == 4 && p.coderate == 1/3
        interleaver.K1 = 860;
        bicm.K1 = 648*p.mu;
        interleaver.K2 = 828;
        bicm.K2 = 624*p.mu;        
    elseif p.Delta_k == 2 && p.coderate == 5/6
        PCCC1.code_param.pun_pattern = [1 0 1 0 1; 1 0 0 1 0; 0 0 0 0 0; 0 0 0 1 0]; % 5/6 coderate
        PCCC2.code_param.pun_pattern = [1 0 1 0 1; 1 0 0 1 0; 0 0 0 0 0; 0 0 0 1 0]; % 5/6 coderate
        bicm.K1 = 624*p.mu;
        bicm.K2 = 576*p.mu;
        interleaver.K1 = 2070;
        interleaver.K2 = 1910;
    elseif p.Delta_k == 3 && p.coderate == 5/6
        PCCC1.code_param.pun_pattern = [1 0 1 0 1; 1 0 0 1 0; 0 0 0 0 0; 0 0 0 1 0]; % 5/6 coderate
        PCCC2.code_param.pun_pattern = [1 0 1 0 1; 1 0 0 1 0; 0 0 0 0 0; 0 0 0 1 0]; % 5/6 coderate
        interleaver.K1 = 2123;
        bicm.K1 = 640*p.mu;
        interleaver.K2 = 2018;
        bicm.K2 = 608*p.mu;
        PCCC2.code_param.tail_pattern(end,end)=0;
    elseif p.Delta_k == 4 && p.coderate == 5/6
        PCCC1.code_param.pun_pattern = [1 0 1 0 1; 1 0 0 1 0; 0 0 0 0 0; 0 0 0 1 0]; % 5/6 coderate
        PCCC2.code_param.pun_pattern = [1 0 1 0 1; 1 0 0 1 0; 0 0 0 0 0; 0 0 0 1 0]; % 5/6 coderate
        bicm.K1 = 648*p.mu;
        bicm.K2 = 624*p.mu;
        interleaver.K1 = 2150;
        interleaver.K2 = 2070;
    else
        msg = 'Undefined interleaver size, Delta_k can be only {1, 2, 4}';
        error(msg)
    end
else 
    defined_KM = 672;
    msg = 'Sorry! cannot initialize turbo coding \nUndefined interleaver size \nInterleaver size is defined only for KM = ';
    error([msg '%d \nHave a good day :-)'],defined_KM)
end

end
