function [ payload, config] = get_payload(config, noblk, mode)

N_data = config.payload.N_data;
Np = config.Np;
nT = config.nT;
nR = config.nR;


if strcmp(mode,'dummy')
    payload = (randn(N_data, nT, noblk)+1i*randn(N_data, nT, noblk))/sqrt(2);    
    
elseif strcmp(mode,'coded')
    p = config.payload.p;    
    
    Kon = p.Kon;
    N_data_on = Kon*p.Mon;
    
    [PCCC, interleaver_size, bicm_size] = find_interleaver_size(p, N_data_on, nT);
        
    dd = cell(noblk,1);
%     dd_short = cell(noblk,1);
    Dd = zeros(p.K,p.M,nT,noblk);
    xd_iT = zeros(N_data,nT,noblk);   
        
    b_blk = cell(noblk,1);
    dd_nT = cell(noblk,1);
    bc_nT_blk = cell(noblk,1);
    code_interleaver = cell(noblk,1);
    bicm_interleaver = cell(noblk,1);
    for blk = 1:noblk
        % Random interleavers
        code_interleaver{blk} = randperm(interleaver_size)-1;
        bicm_interleaver{blk} = randperm(bicm_size)-1;
        
        % Random binary bits
        b_blk{blk} = round(rand([1 interleaver_size]));
        
        % Encoder
        [dd_nT{blk}, bc_nT_blk{blk} ] = do_encode_pccc(b_blk{blk}, PCCC, code_interleaver{blk}, bicm_interleaver{blk});
        
        % Spatial multiplexing
        dd{blk} = reshape(dd_nT{blk},[N_data_on nT]);
                
        for iT = 1:nT
            % Time-Frequency allocation
            Dd(:,:,iT,blk) = do_map(p, dd{blk}(:,iT) );
            
%             Dd(:,:,iT,blk) = 0;
%             Dd(1,1,iT,blk) = 0.99*sqrt(p.Kon);
            
%             Dd(:,1,iT,blk) = 0;
%             Dd(floor(p.K/2),1,iT,blk) = sqrt(p.K);
            
            % GFDM Modulation
            xd_iT(:,iT,blk) = do_modulate(p, Dd(:,:,iT,blk));
            
        end
    end
    config.matrices.Dd_indx = logical(abs(Dd(:,:,:,1)));
    
    [on_bins, on_bins1, off_bins, on_bins_single] = get_on_bins(fft(xd_iT(:,1,1)), N_data, nT,nR);
    config.off_bins = off_bins;
    config.on_bins = on_bins;
    config.on_bins1 = on_bins1;
    config.on_bins_single = on_bins_single;
    
    config.payload.code_interleaver = code_interleaver;
    config.payload.bicm_interleaver = bicm_interleaver;
    config.payload.b_blk = b_blk;
    config.payload.bc_nT_blk = bc_nT_blk;
    config.payload.dd_nT = dd_nT;   
    config.payload.PCCC = PCCC;
    
    payload = xd_iT;
end



end

