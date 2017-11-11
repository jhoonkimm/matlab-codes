function [] = saveFigure(figureFolder, name)

if (~exist('name', 'var'))
  name = figureFolder;
  figureFolder = '';
end

% global shouldSaveFigures
% figureFolder = './figs/';
% if (shouldSaveFigures)
set(gcf, 'PaperPositionMode', 'auto');
%   saveas(gcf, [figureFolder name], 'epsc');
%   print('-depsc2', [figureFolder name]);
%     epswrite([figureFolder name]);
%   print('-dill', [figureFolder name]);
try
  epswrite([figureFolder name]);
catch e
  e
end
end

% end