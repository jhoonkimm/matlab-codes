function [newvar] = interpGaitCycle(var, varargin)

% INTERPGAITCYCLE resample gait cycle to contain N points
%  Resample gait cycle (or other variable) to strecth across N points.
%  If varargin is unspecified, the default number of points (N) is set to 1000.
% 
%  If input VAR contains NaN values, NaN values will also appear in NEWVAR output.
%  NANTOZERO function can be used to convert NaN values to zeros.
%  
%  Example:
%    y=1:100; newy = interpGaitCycle(y); % length(newy) == 1000
%    y=1:100; newy = interpGaitCycle(y,200); % length(newy) == 200
%    y=1:100; newy = interpGaitCycle(y,45); % length(newy) == 45

%  Karl Zelik
%  2/2/10

if nargin==1
    N=1000;
elseif nargin==2
    N=varargin{1};
else
    fprintf('Error: Unexpected varargin in interpGaitCycle.m');
    return;
end


if size(var,2) > size(var,1)
    var = var';
end

for j=1:size(var,2)
    newvar(:,j) = interp1([1:size(var,1)]',var(:,j),linspace(1,length(var),N)); 
end
    

% % Replace NaN values from interpolation with zeros
% nan_indices = find(isnan(newvar)); 
% newvar(nan_indices)=0;