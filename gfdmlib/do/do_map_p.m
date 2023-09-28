% Copyright (c) 2014 TU Dresden
% All rights reserved.
% See accompanying license.txt for details.
%


function D = do_map_p(p, s)
% Map the linear symbol stream to the GFDM data matrix.
%
% Respects the kset and mset values. Fills the matrix sub-symbol
% wise, i.e. the first symbols in the stream belong to the first
% sub-symbol, the next symbols belong to the 2nd and so on.
if (~isfield(p,'mode'))
    p.mode = 'normal';
end

if (p.Delta_k ==1 || p.Mp == 0)
    
    kset = get_kset(p)+1;
    mset = get_mset(p)+1;
    
    if length(kset)==p.K && length(mset) == p.M
        D = reshape(s, [p.K, p.M]);
    else
        Dm = reshape(s, [length(kset), length(mset)]);
        res1 = zeros(p.K, length(mset));
        res1(kset, :) = Dm;
        res = zeros(p.K, p.M);
        res(:,mset) = res1;
        D = res;
    end
    
elseif (strcmp(p.mode,'2blocks') && p.Mp ~= 0)
    D_set = true(p.Kon,p.M);
    %D_set([1:p.k_start end-p.k_end+1:end],:) = false;
    if p.Mp == 2
    D_set(1:p.Delta_k:end,1) = false;
    D_set(1:p.Delta_k:end,end) = false;
    end
    d_set = D_set(:);
    
    s_long = zeros(p.Kon*p.M,1);    
    s_long(d_set) = s;
    
    kset = get_kset(p)+1;
    mset = get_mset(p)+1;
    
    if length(kset)==p.K && length(mset) == p.M
        D = reshape(s_long, [p.K, p.M]);
    else
        Dm = reshape(s_long, [length(kset), length(mset)]);
        res1 = zeros(p.K, length(mset));
        res1(kset, :) = Dm;
        res = zeros(p.K, p.M);
        res(:,mset) = res1;
        D = res;
    end
    
    
elseif strcmp(p.mode,'precode')
    d_set = 1:p.M*p.Kon;
    if p.Mp == 1
    d_set(1:p.Delta_k:p.Kon) = [];
    elseif p.Mp == 2
        d_set([1:p.Delta_k:p.Kon p.Kon+1:p.Delta_k:p.Kon+p.Kon]) = [];        
    end
    
    s_long(d_set) = s;
    
    kset = get_kset(p)+1;
    mset = get_mset(p)+1;
    
    if length(kset)==p.K && length(mset) == p.M
        D = reshape(s_long, [p.K, p.M]);
    else
        Dm = reshape(s_long, [length(kset), length(mset)]);
        res1 = zeros(p.K, length(mset));
        res1(kset, :) = Dm;
        res = zeros(p.K, p.M);
        res(:,mset) = res1;
        D = res;
    end
else    
    d_set = 1:p.M*p.K;
    if p.Mp == 1
    d_set(1:p.Delta_k:p.Kon) = [];
    elseif p.Mp == 2
        d_set([1:p.Delta_k:p.K end-p.K+1:p.Delta_k:end]) = [];        
%         d_set([p.Kon+1:p.Delta_k:p.Kon+p.Kon end-p.Kon+1:p.Delta_k:end]) = [];        
    end
    
    s_long(d_set) = s;
    
    kset = get_kset(p)+1;
    mset = get_mset(p)+1;
    
    if length(kset)==p.K && length(mset) == p.M
        D = reshape(s_long, [p.K, p.M]);
    else
        Dm = reshape(s_long, [length(kset), length(mset)]);
        res1 = zeros(p.K, length(mset));
        res1(kset, :) = Dm;
        res = zeros(p.K, p.M);
        res(:,mset) = res1;
        D = res;
    end
end
