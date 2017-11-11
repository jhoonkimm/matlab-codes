function out = vecnorm(data,varargin)

% VECNORM computes vector norm along 1 dimension of a matrix
%  (i.e., computes magnitude of a vector)
% 
% OUT = VECNORM(data,dim,ord)
%   data: vector/matrix of data
%   dim: vectors are defined across this dimenion
%        (e.g., dim=2 means each row is defined as a vector)
%   ord: vector norm order; 2nd order vetor norm calculates
%        magnitude of a vector. 
% 
% Take the norm of a matrix, treating groups along dimension "dim" 
% as vectors to be normed. Order of norm is "ord". 
% Default is a 2nd order norm along groups in the smallest 
% non-singleton dimension. 
% 
% Examples: 
% x = [1 2; ...
%      3 4; ...
%      5 6];
% 
% vecnorm(x,2,2) 
% ans =  [5; 25; 61].^(1/2)
%     2.2361
%     5.0000
%     7.8102
%  
% vecnorm(x,1,2)
% ans = [35 56].^(1/2)
%     5.9161    7.4833  

% Peter Adamczyk


if nargin == 1
    dim = size(data);
    dim = min(find(dim~=1)) ; % 
    ord = 2;
elseif nargin == 2
    dim = varargin{1} ;
    ord = 2 ;
elseif nargin == 3
    dim = varargin{1};
    ord = varargin{2};
else
    error('01-nargin','Too Many Inputs')
end

if ord==inf
    out = max(data,[],dim);
elseif ord == -inf
    out = min(data,[],dim);
% elseif ord == 2
%     out = sqrt(sum(abs(data).^ord,dim)) ;
else
    out = sum(abs(data).^ord,dim).^(1/ord) ;
end

    
    
    