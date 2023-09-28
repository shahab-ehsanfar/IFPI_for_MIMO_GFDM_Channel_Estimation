% Copyright (c) 2014 TU Dresden
% All rights reserved.
% See accompanying license.txt for details.
%


function s = do_unmap_p(p, D)
% Convert data matrix to linear symbol stream
%
% Respects the values of kset and mset. The first column of D
% constitutes the first symbols of s, using only these rows that
% exist in kset.

if p.Delta_k == 1
    
    kset = get_kset(p)+1;
    mset = get_mset(p)+1;
    
    if length(kset) == p.K && length(mset) == p.M
        s = reshape(D, numel(D), 1);
    else
        Dm = D(kset,mset);
        s = reshape(Dm, numel(Dm), 1);
    end
    
else
   
    d_set = 1:p.M*p.Kon;
    if p.Mp == 1
    d_set(1:p.Delta_k:p.Kon) = [];
    elseif p.Mp == 2
        d_set([1:p.Delta_k:p.Kon p.Kon+1:p.Delta_k:p.Kon+p.Kon]) = [];        
    end
    
    kset = get_kset(p)+1;
    mset = get_mset(p)+1;
    
    if length(kset) == p.K && length(mset) == p.M
        s_long = reshape(D, numel(D), 1);
    else
        Dm = D(kset,mset);
        s_long = reshape(Dm, numel(Dm), 1);
    end
    
    s = s_long(d_set);   
end

end