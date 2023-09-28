function [ coded_sym,codeword ] = do_encode_pccc(data, struct, code_interleaver, bicm_interleaver)

code_param = struct.code_param;
sim_param = struct.sim_param;

codeword = TurboEncode( data, code_interleaver, code_param.pun_pattern, code_param.tail_pattern, sim_param.g1, sim_param.nsc_flag1, sim_param.g2, sim_param.nsc_flag2 );

codeword = Interleave( codeword, bicm_interleaver );

s1 =  Modulate( codeword, code_param.S_matrix );

coded_sym = s1.';

end
