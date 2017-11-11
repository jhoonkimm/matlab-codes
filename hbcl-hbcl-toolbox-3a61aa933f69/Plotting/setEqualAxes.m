function [] = setEqualAxes(figs)

% SETEQUALAXES rescale two or more figures to have identical axes
% 
%  Warning: Rescaling axes will remove legends, so use SETEQUALAXES
%  before adding legends
% 
%  Examples:
%   figure(1); plot(..); figure(2); plot(..);  setEqualAxes([gca(1) gca(2)])
%   
%   figure; G=subplot(331); plot(..); H=subplot(332); plot(..); setEqualAxes([G H])

% Karl Zelik
% updated 2/11/09

array = [];
for i=1:length(figs)
    n = figs(i);
    oldaxes = axis(n);
    array = [array; oldaxes];
end

newaxes = [min(array(:,1)) max(array(:,2)) min(array(:,3)) max(array(:,4))];

for i=1:length(figs)
    n=figs(i);
    axes(n); axis(newaxes);
end
