function [customPhasespaceLog] = processCustomPhasespaceLog(directoryName, varargin)
% Imports treadmill custom log marker data
% 
% Original by John Rebula
% Last updated by Xiao-Yu Fu 6/21/2012


createCachedFile = 1;
forceReprocess = 0;
verbose = 0;

%% Process input arguments
opt_argin = varargin;
while length(opt_argin) >= 2,
  opt = opt_argin{1};
  val = opt_argin{2};
  opt_argin = opt_argin(3:end);
  switch opt
    case 'verbose'
      verbose = val;
    case 'forceReprocess'
      forceReprocess = val;
    case 'createCachedFile'
      createCachedFile = val;
      %     otherwise
      %       error('\nError: unrecognized parameter ''%s''\n',opt)
  end
end

%% Create directory paths
dotLocations = findstr(directoryName, '.');
outputFilename = [directoryName(1:(dotLocations(end) - 1)) '.mat'];

tic
if (forceReprocess || ~exist(outputFilename, 'file'))
  %% Open file, get markers and frame position 
  fid = fopen(directoryName);
  data = textscan(fid, '%s %f %f %f');
  fclose(fid);
  
  datacols=[data{2:4}];
  
  allMarkers = unique(data{1}(strncmp(data{1},'M',1))); % find all markers
  [dummy,sortOrder] = sort(str2double(strrep(allMarkers,'M',''))); % sort them prettily so John doesn't complain %% thanks! (John)
  allMarkers = allMarkers(sortOrder);
  frameIndex = ~strncmp(data{1},'M',1);
  frameBreaks = find(frameIndex);
  numFrames = length(frameBreaks);
  frameBreaks = [frameBreaks; length(data{1})+1];
  
  %% Time code, frame numbers
  firstTimingCol = str2double(data{1}(frameIndex));
  timingAndTriggerData = [firstTimingCol datacols(frameIndex,:)]+48; % weird offset
  customPhasespaceLog.frameNumber = (timingAndTriggerData*[1 256 256^2 0]')'; % transpose to maintain old order
  
  %% Marker data
  for m=1:length(allMarkers)
    markerName = allMarkers{m};
    markerData = nan(numFrames,3);
    
    markerIndex = strcmp(data{1},allMarkers{m});
    markerIndexPositions = find(markerIndex);
    
    [dummy,framePositions] = histc(markerIndexPositions,frameBreaks);
    markerData(framePositions,:) = datacols(markerIndex,:);
    customPhasespaceLog.(markerName) = markerData;
  end
  
  if (verbose)
    subplot(3,1,1);
    plot(timingAndTriggerData)
    title('timing bytes and trigger byte');
    subplot(3,1,2);
    plot(customPhasespaceLog.frameNumber)
    title('frame number from timing bytes');
    subplot(3,1,3);
    plot(diff(customPhasespaceLog.frameNumber))
    title('diff(frame number from timing bytes)');
  end
  
  if (createCachedFile)
    save(outputFilename, 'customPhasespaceLog');
  end
  
else
  load(outputFilename);
end


end
