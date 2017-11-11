function [] = addFigureFrameToAVI(aviWriter, varargin)
%%
shouldAntiAlias = 0;
figureHandleToSave = [];
for i = 1 : 2 : length(varargin)
  option = varargin{i};
  value = varargin{i + 1};
  switch option
    case 'shouldAntiAlias'
      shouldAntiAlias = value;
    case 'figureHandleToSave'
      figureHandleToSave = value;
  end
end

%%
if (~isempty(figureHandleToSave))
  figure(figureHandleToSave);
end

%%
if (shouldAntiAlias)
  antiAliasedFig = myaa();
  drawnow;
  pause(0.01);
  thisFrame = getframe(antiAliasedFig);
  close(antiAliasedFig);
else
  thisFrame = getframe(gcf);
end

writeVideo(aviWriter, thisFrame);

end


