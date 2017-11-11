function [c3dTrial] = renameMarkersInC3DData(c3dTrial, markerRenaming)

if (~isfield(c3dTrial, 'markerNames'))
  c3dTrial.markerNames = fieldnames(c3dTrial);
end

markerNamesOriginal = c3dTrial.markerNames;
c3dTrial = rmfield(c3dTrial, 'markerNames');

for i = 1:length(markerNamesOriginal)
  renamed = 0;
  for j = 1:length(markerRenaming)
    if (strcmp(markerRenaming{j}{1}, markerNamesOriginal{i}))
      
      fprintf('renaming marker %s to %s\n', markerRenaming{j}{1}, markerRenaming{j}{2});
      
      [c3dTrial2.(markerRenaming{j}{2})] = c3dTrial.(markerRenaming{j}{1});
      %         if (isfield(trackedC3dTrial2, markerRenaming{j}{1}))
      %           trackedC3dTrial2 = rmfield(trackedC3dTrial2, markerRenaming{j}{1});
      %         end
      renamed = 1;
      break;
    end
  end
  if (~renamed)
    [c3dTrial2.(markerNamesOriginal{i})] = c3dTrial.(markerNamesOriginal{i});
  end
end
c3dTrial = c3dTrial2;
c3dTrial.markerNames = fieldnames(c3dTrial);

end

