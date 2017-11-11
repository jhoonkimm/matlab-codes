function [rigidBodySegments, segmentNames] = parseRigidBodiesFromVisual3DModelTemplateFile(templateFilename)

segmentNames = {};
% templateFile = [directory 'modelTemplate.mdh'];
fid = fopen(templateFilename);
block_size = 1;
formatString = '%s %s %s';

index = 1;
while ~feof(fid)
  segarray = textscan(fid, formatString, block_size);
  
  strings = segarray;
  
  if (length(strings) == 3)
    if (isempty(strings{1}) ||  ...
        isempty(strings{2}) || ...
        isempty(strings{3}))
      continue;
    end
    
    if (strcmp(strings{1}{:}, '!') && strcmp(strings{2}{:}, 'Segment') && ~strcmp(strings{3}{:}, 'Info'))
      segmentNames{index} = strings{3}{:};
      index = index + 1;
    end
  end
end
% fclose(fid);

fid2 = fopen(templateFilename);
% frewind(fid)

formatString = '%s';
index = 1;
while ~feof(fid2)
  segarray = textscan(fid2, formatString, block_size);
  strings = segarray;
  
  if (length(strings) == 1)
    if (isempty(strings{1}))
      continue;
    end
    
    [ind] = findstr(strings{1}{:}, '/TRACKING_NAMES=');
    
    if (~isempty(ind))
      %       [ind] = findstr(strings{1}{:}, segmentString);
      markerNames = strings{1}{:};
      
      markerNames(1:16) = [];
      markerNameList = regexp(markerNames, '+', 'split');
      rigidBodySegments.(segmentNames{index}) = markerNameList;
      
      index = index + 1;
    end
  end
end

fclose(fid);


end

