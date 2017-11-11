function [varargout] = breakout(varargin)

% BREAKOUT break large matrix into individual column or row arrays
%  This function breaks out data from a big matrix into all its components,
%  maintaining the integrity of each subunit appearing as a group with
%  constant dir'th index (dir for "direction")
% 
%  VARARGIN = {mat, dir, inds1, inds2, ...}
%  Function reuqires at least 2 input arguments (mat & dir)
% 
% Examples:  
%  [a b] = breakout([1 2; 3 4],2)
%    result: a = [1; 3]; b = [2; 4]
% 
%  [a b] = breakout([1 2; 3 4],1)
%    result: a = [1 2]; b = [3 4]

%  Peter Adamczyk


mat = varargin{1};
dir = varargin{2};
lastgot = 0;
ndims = length(size(mat));

for i = 1:nargout
    S{i}.type = '()';
    for j = 1:ndims
        S{i}.subs{j} = ':' ;
        if j == dir
            try
                S{i}.subs{j} = varargin{i+2};
                lastgot = max(varargin{i+2});
            catch
                S{i}.subs{j} = lastgot+ [1 : ceil((size(mat,dir)-lastgot)/(nargout-i+1))]  ;
                lastgot = max(S{i}.subs{j});
            end
        end
    end
end

for i = 1:nargout
    varargout{i} = squeeze(subsref(mat,S{i}));
end
