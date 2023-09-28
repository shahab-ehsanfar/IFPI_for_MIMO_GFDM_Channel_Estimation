% Copyright (c) 2014 TU Dresden
% All rights reserved.
% See accompanying license.txt for details.
%


function P = get_power_delay_profile(numPaths,nT,nR)

P = zeros(numPaths,nT,nR);
for iR = 1:nR
    for iT = 1:nT
        AveGaindb = linspace(0,-20, numPaths);
        a = 10.^(AveGaindb/10);
        P(:,iT,iR) =  a./(sum(a)); % a./sqrt(sum(a.^2));
    end
end

end
