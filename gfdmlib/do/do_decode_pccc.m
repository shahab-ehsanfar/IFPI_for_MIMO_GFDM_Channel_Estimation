function [ detected_data, bit_likelihood ] = do_decode_pccc(sym_hat, data, EsNo, a, struct, code_interleaver, bicm_interleaver)

code_param = struct.code_param;
sim_param = struct.sim_param;
turbo_iterations = 8;

%a = ones(1,code_param.symbols_per_frame);
%EsNo = 1000;
input_somap_c = zeros(1, code_param.code_bits_per_frame );

symbol_likelihood = Demod2D( sym_hat.', code_param.S_matrix, EsNo, a );
%sim_param.demod_type = 4;

if (isfield(struct,'hard') && struct.hard == true)    
    bit_likelihood = Somap( symbol_likelihood, sim_param.demod_type, input_somap_c );
    detected_data = (sign(bit_likelihood)+1)/2;
else
    for bicm_iter = 1:8
        % demodulate
        bit_likelihood = Somap( symbol_likelihood, sim_param.demod_type, input_somap_c );
        input_decoder_c = bit_likelihood(1:code_param.code_bits_per_frame);
        
        % deinterleave (bicm)
        input_decoder_c = Deinterleave( input_decoder_c, bicm_interleaver);
        
        % decode
        [detected_data, turbo_errors, output_decoder_c] = TurboDecode( input_decoder_c, data, turbo_iterations, sim_param.decoder_type, code_interleaver, code_param.pun_pattern, code_param.tail_pattern, sim_param.g1, sim_param.nsc_flag1, sim_param.g2, sim_param.nsc_flag2);
        
        
    end    
end

end
