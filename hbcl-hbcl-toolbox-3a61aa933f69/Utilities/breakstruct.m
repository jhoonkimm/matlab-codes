function varargout = breakstruct(varargin)

% BREAKSTRUCT grab structure fields and store as workspace variables  
%
% breakstruct(data,flds,dim) breaks out fields of the input array "data" 
% into their respective arrays as named in "flds", storing the newly 
% created variables In The Calling Workspace.  
%   To create an array, the field values are concatenated according to
%   a single dimension, with this priority:
%         (1) dim (specified input) if available
%         (2) first non-singleton dimension of the input array (if field value is scalar)
%         (3) first _singleton_ dimension of the field value array (if vector)
%         (4) a new dimension if all other dimensions are non-singleton
% 
% out = breakstruct(data,flds,dim)  breaks out the specified field(s) and
% returns them to the assigned variables. 
% 
% "data" is a structure array (array of structures) 
% that includes fields "flds", and possibly other fields as well
% "flds" is a cell array of strings, such as is returned by fieldnames,
% or created as in:
%   flds = {'t','x','y'};
% dim is optional, specifying the dimension along which each field's data
% will be concatenated
% 
% if "flds" is not specified, ALL fields of the input structure are used
% 

% if nargin == 1
%     data = varargin{1} ;
%     flds = fieldnames(data) ;
%     dimspec = [];
% elseif nargin == 2
%     data = varargin{1} ;
%     flds = varargin{2} ;
%     dimspec = [];
% elseif nargin == 3
%     data = varargin{1} ;
%     flds = varargin{2} ;
%     dimspec = varargin{3} ;
% else
%     errmsg('Error::BREAKSTRUCT::inputs','Too many inputs to BREAKSTRUCT') ;
% end

% default is no dimspec
flds = {};
dimspec = [];

% now overwrite if they are in a different order or if "flds" or "dim" is
% supplied
for i = 1:nargin
    if isstruct(varargin{i})
        data = varargin{i} ;
        continue
    elseif iscell(varargin{i})
        flds = varargin{i} ;
        continue
    elseif isnumeric(varargin{i})
        dimspec = varargin{i} ;
        continue
    else
        error(['Error::BREAKSTRUCT::Unacceptable Input Type, input ' int2str(i)],'Unknown Input Type') ;
    end
end
if isempty(flds)
    flds = fieldnames(data);
end

arrdim = find(size(data)>1, 1, 'first') ; % first non-singleton dimension of "data"
if isempty(arrdim) % scalar?!
    arrdim = 1;
end

numflds = length(flds) ;
for i = 1:numflds
    
    if ~isempty(dimspec) % if dimspec exists, then it overrides all other dimension choices
        dim = dimspec ;
    elseif isscalar(data(1).(flds{i}))
% 2007-09-06 PGA        dim = arrdim ; % if this field in each structure is a scalar, concatenate into the same shape as the original array
        dim = size(data);
    elseif any( size(data(1).(flds{i})) == 1 )
        dim = min(find(size(data(1).(flds{i}))==1)) ; % first non-singleton dimension of the FIELD value array
    else % if we get here, then all dimensions are non-singleton and we must concatenate in an additional dimension. 
        dim = length(size(data(1).(flds{i})))+1 ;
    end
    
    if ~isscalar(dim)  % if "dim" has more than one element, it is the size of the intended result.
        temp = arrayfun(@(x) x.(flds{i}), data);
    else
        temp = cat(dim,data.(flds{i})) ;
    end
    
    if nargout == 0
        assignin('caller',flds{i},temp);
    else 
        outputs{i} = temp ;
    end
end

if nargout == 0
    return
elseif nargout==1
    varargout{1} = cat(find(~ismember([1:10], [arrdim]),1,'first'), outputs{:}) ;
else
    varargout = outputs ;
end

end

