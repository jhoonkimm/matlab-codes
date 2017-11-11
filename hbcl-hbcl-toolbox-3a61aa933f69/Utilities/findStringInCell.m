function [inds] = findStringInCell(cellin,strin,varargin)

% FINDSTRINGINCELL find a string within a cell
%   
%  Example:
%   temp = {'all','your','base','are','belong','to','us'};
%   findStringInCell(temp,'base') % returns 3
%   findStringInCell(temp,'touchdown') % returns []
%   findStringInCell(temp,'o') % returns []

%   Karl Zelik 11/19/09

% Default varaible
wholeWord = 1; % find exact match to entire input string?

% Optional input arguments
opt_argin = varargin;
while length(opt_argin) >= 2,
  opt = opt_argin{1};
  val = opt_argin{2};
  opt_argin = opt_argin(3:end);
  switch opt
    case 'wholeWord'
      wholeWord = val;
    otherwise
      error('\nError: incorrect varargin\n')
  end
end

inds = [];
if ~wholeWord
    for i = 1:length(cellin)
        if ~isempty(findstr(cellin{i},strin))
            inds = [inds; i];
        end
    end
else
    for i = 1:length(cellin)
        if strcmp(cellin{i},strin)
            inds = [inds; i];
        end
    end
end