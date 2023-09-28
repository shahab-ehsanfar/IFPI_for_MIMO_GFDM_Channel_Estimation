function [ y_out ] = do_emulate_cp_circUW(config,h_hat,y,blk,mode,xd)
% This function emulates the CP at the receiver side
% In other words, the channel h_hat looks circulant to the signal of interest
%
% TU Dresden
% Shahab

nT = config.nT;
nR = config.nR;
xp_iT = config.preamble.xp_iT;
N_uw = config.Np;
N_data = config.payload.N_data;


if strcmp(mode,'for_data')
    %% CP Emulation for data 
    numPaths = config.numPaths; 
    
    % Copy the UW part and add it to the beginning of the block
    y_d0 = zeros(N_data,nR);
    for iR = 1:nR
        y_d0(:,iR) = y(N_uw+1:N_uw+N_data,iR,blk) + padarray(y(N_uw+N_data+1:N_uw+N_data+numPaths,iR,blk),N_data-numPaths,'post');
    end
    
    % Calculate the signal from the known terms
    xp_interf1 = zeros(N_uw+numPaths,nR);
    xp_interf2 = zeros(N_uw,nR);
    xp_interf = zeros(numPaths,nR);
    if (mod(blk,2))
        for iR = 1:nR
            for iT = 1:nT
                xp_interf1(:,iR) = xp_interf1(:,iR) + filter(h_hat(:,iT,iR),1,[xp_iT(N_uw+1:end,iT); zeros(numPaths,1)]); % interference of UW to UW-Payload
                xp_interf2(:,iR) = xp_interf2(:,iR) + filter(h_hat(:,iT,iR),1, xp_iT(N_uw+1:end, mod((iT-1)+2,nT)+1 ) ); % interference from copying UW in Payload-UW to beginning
                xp_interf(:,iR) = xp_interf1(N_uw+1:end,iR) + xp_interf2(1:numPaths,iR);
            end
        end
    else
        for iR = 1:nR
            for iT = 1:nT
                xp_interf1(:,iR) = xp_interf1(:,iR) + filter(h_hat(:,iT,iR),1,[xp_iT(N_uw+1:end, mod((iT-1)+2,nT)+1 ); zeros(numPaths,1)]); % interference of UW to UW-Payload
                xp_interf2(:,iR) = xp_interf2(:,iR) + filter(h_hat(:,iT,iR),1, xp_iT(N_uw+1:end, iT ) ); % interference from copying UW in Payload-UW to beginning
                xp_interf(:,iR) = xp_interf1(N_uw+1:end,iR) + xp_interf2(1:numPaths,iR);
            end
        end
    end
    
    % Remove the uw-interference and extract the circulant-channeled data signal
    yd = zeros(N_data,nR);
    for iR = 1:nR
        yd(:,iR) = y_d0(:,iR) - padarray(xp_interf(1:numPaths,iR),N_data - numPaths,'post');
    end    
    y_out = yd;

% elseif strcmp(mode,'for_uw')
%     %% CP Emulation for the Unique Word (Pilots) part
% 
%     numPaths = config.numPaths;
%     
%     % Copy the beginning of data part and add it to the beginning of the UW part
%     y_uw0 = zeros(N_uw,nR);
%     for iR = 1:nR
%         y_uw0(:,iR) = y(1:N_uw,iR,blk) + y(N_data+1:end,iR,blk); 
%     end    
%        
%     % Calculate the ISI leakage from data to UW part
%     x_leakage = zeros(numPaths-1+N_uw,nT);
%     x_padded = zeros(2*(numPaths-1)+N_uw,nT);
%     for iT = 1:nT
%         x_leakage(:,iT) = [xd(end-numPaths+2:end,iT); xd(1:N_uw,iT)];
%         x_padded(:,iT) = padarray(x_leakage(:,iT),numPaths-1,'pre');    
%     end
%     x_filtered = zeros(2*(numPaths-1)+N_uw,nR);
%     interference = zeros(N_uw,nR);
%     for iR = 1:nR
%        for iT = 1:nT
%            x_filtered(:,iR) = x_filtered(:,iR) + filter(h_hat(:,iT,iR),1,x_padded(:,iT));              
%        end
%        interference(:,iR) = x_filtered(end-N_uw+1:end,iR);
%     end   
%         
%     % Remove the interference and extract the circulant-channeled UW signal
%     y_uw = zeros(N_uw,nR);
%     for iR = 1:nR
%        y_uw(:,iR) = y_uw0(:,iR) - interference(:,iR);
%     end    
%     y_out = y_uw;
        
    
end





end

