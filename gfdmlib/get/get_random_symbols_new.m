% Copyright (c) 2014 TU Dresden
% All rights reserved.
% See accompanying license.txt for details.
%


function s = get_random_symbols_new(p)
% Return linear stream of random symbols

if (p.Delta_k ==1)
    
    num = length(get_mset(p)) * length(get_kset(p));
    
    high = 2 ^ p.mu;
    
    s = randi(high, [num, 1])-1;

else    
    p0 = p; p0.M= p.Md+ (p.Delta_k-1)*p.M; 
    p0.Kset = 1:p.Kon/p.Delta_k;
        
    num = length(get_mset(p0)) * length(get_kset(p0));
    
    high = 2 ^ p.mu;
    
    s = randi(high, [num, 1])-1;   
    
end


