function [ struct ] = structurize(varargin)

for i = 1:nargin   
    eval(['struct.' inputname(i) '= varargin{i};'])        
end

end