function buildStructure(name, varargin)

% BUILDSTRUCTURE create or append to structure from workspace variables
%  Save workspace variables into structure
%  
%  BUILDSTRUCTURE(NAME, VARARGIN)
%   Argument NAME should be input as STRING.
%
%  VARARGIN Options:
%   vars: specify which workspace variables to save into structure;
%         argument should be a cell of strings (e.g., {'x','y','z'})
%         default is to save all workspace variables
%   overwrite: overwrite (1) or append (0) to structure;
%         default is to append to structure with NAME if it exists
%
%  Example:
%    clear; x=2; y=5; z=10;
%    buildStructure('Q');
%    result: Q=1x1 struct with fields x,y,z
%
%    buildStructure('Q');
%    result: now Q=1x2 struct with fields x,y,z
%
%    buildStructure('Q','vars',{'y','z'},'overwrite',1)
%    result: now Q=1x1 struct with fields y,z only

% Karl Zelik 11/22/09
% based on code (stackdata.m, buildstruct.m) from Peter Adamczyk

% Default variables
vars = []; % if empty, all worksapce variables are saved
overwrite = 0; % by default, append to stucture rather than overwrite it

% Optional input arguments
opt_argin = varargin;
while length(opt_argin) >= 2,
  opt = opt_argin{1};
  val = opt_argin{2};
  opt_argin = opt_argin(3:end);
  switch opt
    case 'vars'
      vars = val;
    case 'overwrite'
      overwrite= val;
    otherwise
      error('\nError: incorrect varargin\n')
  end
end

% If overwrite structure
if overwrite
   evalin('base',['clear ' name]) 
end

% If vars is empty, save all variables in workspace 
if isempty(vars)
    vars = evalin('caller','who');
end

try
    sz = evalin('caller',['length(' name ')']) ;
catch
    sz = 0;
end

ind = ['(' int2str(sz+1) ')'];

for ii = 1:length(vars)
    if ~strcmp(vars{ii},name)
        evalin('caller',[name ind '.' vars{ii} ' = ' vars{ii} ';']);
    end
end