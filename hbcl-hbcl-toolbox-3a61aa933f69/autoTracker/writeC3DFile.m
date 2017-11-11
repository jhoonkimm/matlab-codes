function [] = writeC3DFile(outputFilename, allMarkerData, varargin)
%%


analogData = [];
analogToVideoRatio = 0; %analogToMarkerRatio;
videoFrameRate = 480;

for i = 1 : 2 : length(varargin)
  option = varargin{i};
  value = varargin{i + 1};
  switch option
    case 'analogData'
      analogData = value;
    case 'analogToVideoRatio'
      analogToVideoRatio = value;
    case 'videoFrameRate'
      videoFrameRate = value;
  end
end

numberOfAnalogChannels = 0; %length(analogChannelData);
if ~isempty(analogData)
  analogChannelNames = fieldnames(analogData);
  numberOfAnalogChannels = length(analogChannelNames);
  for i = 1 : length(analogChannelNames)
    analogChannelData{i} = analogData.(analogChannelNames{i});
  end
end



%% don't ever ever name variables like this, name them something sane.
%   vfrt     = video frame rate, must be > 0, default 120
%   nfram    = number of frames, default 500
%   nmkr     = number of markers, not < 0, default 20
%   avr      = analog to video ratio, default 13
%   achn     = number of analog channels, not < 0, default 28
%   ftype    = Intel(1), DEC(2), SGI(3), default 1
%   dtype    = Integer(1), Floating Point(2), default 2
%   pscal    = point data scaling factor, must be > 0, default 0.1

c3dServer = actxserver('C3DServer.C3D');
name = outputFilename; %'output.c3d';

markerNames = fieldnames(allMarkerData);
numberOfFrames = length(allMarkerData.(markerNames{1}));
filetype = 3; %1; %
dataType = 2; %
pointDataScaleFactor = 1.0;
numberOfMarkers = length(markerNames);
errorOnFindingNaNs = 0;

createc3d(c3dServer, name, videoFrameRate, numberOfFrames, numberOfMarkers, ...
    analogToVideoRatio, numberOfAnalogChannels, filetype, dataType, pointDataScaleFactor);


%% write analogs:

if ~isempty(analogData)
  nIndex = c3dServer.GetParameterIndex('ANALOG', 'OFFSET');
  c3dServer.SetParameterType(nIndex, 4);
  analogIndex = 0;
  for analogNumber = 1:numberOfAnalogChannels
    c3dServer.SetParameterValue(nIndex, analogIndex, single(0.0));
    analogIndex = analogIndex + 1;
  end
  
  
  for channelNumber = 0:(numberOfAnalogChannels - 1)
    %   fprintf('writing channel %g\n', channelNumber);
    
    dataToWrite = zeros(1, numberOfFrames * 2);
    dataToWrite(1:length(single(analogChannelData{channelNumber + 1}))) = ...
      single(analogChannelData{channelNumber + 1});
    
    %   subFrameNums = 1:(numberOfFrames * analogToVideoRatio);
    c3dServer.SetAnalogDataEx(channelNumber, 1, single(dataToWrite + 0)); %
    % sin(subFrameNums/10)));
    
    %   frameNums = 1:(numberOfFrames);
    %   for frameNum = frameNums
    %     for subFrameNum = 1:analogToVideoRatio
    %       c3dServer.SetAnalogData(channelNumber, frameNum, subFrameNum, sin(frameNum/10));
    %     end
    %   end
  end
end

%% specify units of marker data:
parameterName = 'UNITS';
description = '';
parameterGroup = 'POINT';
lockParam = '0';
dataType = -1; %-1 : Character, 1 : Byte, 2 : Integer, 4 : Floating point
numberOfDataDimensions = 2; %needs to be 2 for a string, since first dimension holds the length of each string (!?)
dataDimensions = [4, numberOfMarkers]; % [lengthOfEachString, numberOfStrings]
data = {};
for markerNumber = 0:(numberOfMarkers - 1)
  data{length(data) + 1} = 'mm';
end

addParametersToC3D = 1; %0; %
if (addParametersToC3D)
  c3dServer.AddParameter(parameterName, description, ...
    parameterGroup, lockParam, dataType, numberOfDataDimensions, ...
    dataDimensions, data); %[]); %
end

%% add extra parameters

% for i = 1:length(parameterDataToWrite)
%   parameterGroup = 'EXTRA';
%   if (i == 1)
%     c3dServer.AddGroup (0, parameterGroup, 'extra parameters added by matlab exporter', '0');
%   end
%   parameterName = parameterDataToWrite{i}{1}; %'ANGLE';
%   description = '';
%   lockParam = '0';
%   dataType = int16(4); %2); %4; %-1; %-1 : Character, 1 : Byte, 2 : Integer, 4 : Floating point
%   numberOfDataDimensions = int16(2); %int16(2); %2; %2; %needs to be 2 for a string, since first dimension holds the length of each string
%   dataDimensions = int16([1 1]); %1; %{uint8(1)}; %[1 1]; %[1, numberOfMarkers]; % [lengthOfEachString, numberOfStrings]
%   data = single([parameterDataToWrite{i}{2} 0]);
%   
%   addParametersToC3D = 1; %0; %
%   if (addParametersToC3D)
%     c3dServer.AddParameter(parameterName, description, ...
%       parameterGroup, lockParam, dataType, numberOfDataDimensions, ...
%       dataDimensions, data);
%   end
% end



%% name markers
nIndex = c3dServer.GetParameterIndex('POINT', 'LABELS');
% nItems = c3dServer.GetParameterLength(nIndex);
% target_name = c3dServer.GetParameterValue(nIndex, 0);
markerIndex = 0;
for markerNumber = 1:numberOfMarkers
  markerName = markerNames{markerNumber};
  if (strcmp(markerName, 'frameNumber'))
    continue;
  end
  c3dServer.SetParameterValue(nIndex, markerIndex, markerName);
  markerIndex = markerIndex + 1;
end

%% name analog signals
if ~isempty(analogData)
  nIndex = c3dServer.GetParameterIndex('ANALOG', 'LABELS');
  for analogChannelNumber = 0:(length(analogChannelNames) - 1)
    analogChannelName = analogChannelNames{analogChannelNumber + 1}; %(((end-5):end)); %
    %   fprintf('setting analog channel name: %s\n', analogChannelName);
    c3dServer.SetParameterValue(nIndex, analogChannelNumber, analogChannelName);
  end
end

%% write marker data
markerIndex = 0;
for markerNumber = 1 : numberOfMarkers
  markerName = markerNames{markerNumber};
  if (strcmp(markerName, 'frameNumber'))
    continue;
  end
  markerData = allMarkerData.(markerName);
  
  if (sum(sum(isnan(markerData))))
    if (errorOnFindingNaNs)
      error('shoulnd''t output nan data from filled in c3d!');
    end
  end
  
  for channel = 0:3
    frameNums = 1:numberOfFrames;
    if (channel == 3) % residual channel... very important that this be non-negative
      if (sum(isnan(markerData(:, 1)')))
        if (errorOnFindingNaNs)
          error('about to write a nan bit of data from marker %s\n', markerName);
        end
      end
      dataToWrite = ones(size(frameNums));
      dataToWrite(isnan(markerData(:, 1)')) = -1;
      c3dServer.SetPointDataEx(markerIndex, channel, 1, single(dataToWrite));
    else % x, y, z channels
      dataToWrite = markerData(:, channel + 1)';
      if (sum(isnan(dataToWrite)))
        if (errorOnFindingNaNs)
          error('about to write a nan bit of data from marker %s\n', markerName);
        end
      end
      
      channelToWrite = channel;
      % We don't want to do this transformation here
      %       % we want to do a transformation:
      %       % y->z
      %       % -z->y
      %       % x->x
      %       if (channelToWrite == 0)
      %         channelToWrite = 0;
      %       elseif (channelToWrite == 1)
      %         channelToWrite = 2;
      %       elseif (channelToWrite == 2)
      %         channelToWrite = 1;
      %         dataToWrite = -1.0 * dataToWrite;
      %       end
      
      %       if(fillMarkerGapsWithSplines)
      %         goodIndeces = ~isnan(dataToWrite);
      %         allIndeces = 1:length(dataToWrite);
      %         dataToWrite = spline(allIndeces(goodIndeces), dataToWrite(goodIndeces), allIndeces);
      %       end
      
      %       sum(isnan(dataToWrite))
      dataToWrite(isnan(dataToWrite)) = 0;
      
      % This may need to be put back in.
      %       if(c3dUnitsAreInMeters)
      %         markerDataOutputScale = 0.001;
      %       else
      markerDataOutputScale = 1;
      %       end
      
      retVal = c3dServer.SetPointDataEx(markerIndex, channelToWrite, 1, single(dataToWrite * markerDataOutputScale));
      %       retVal
    end
    
    %     for frameNum = frameNums
    %       if (channel == 3) % residual channel... very important that this be non-negative
    %         c3dServer.SetPointData(markerNumber, channel, frameNum, 0.0001);
    %       else % x,y,z channels
    %         c3dServer.SetPointData(markerNumber, channel, frameNum, sin(frameNum/10));
    %       end
    %     end
  end
  markerIndex = markerIndex + 1;
end


%%
% numberOfFrames = numberOfFrames - 2; %!?
parameterGroup = 'POINT';
parameterName = 'LONG_FRAMES';
description = '';
lockParam = '0';
dataType = int16(4); %2); %4; %-1; %-1 : Character, 1 : Byte, 2 : Integer, 4 : Floating point
numberOfDataDimensions = int16(2); %int16(1); %2; %2; %needs to be 2 for a string, since first dimension holds the length of each string
dataDimensions = int16([1 1]); %1; %{uint8(1)}; %[1 1]; %[1, numberOfMarkers]; % [lengthOfEachString, numberOfStrings]
data = single([numberOfFrames 0]); %int16([mod(numberOfFrames, (2^16 - 1)), floor(numberOfFrames / (2^16 - 1))]);

c3dServer.AddParameter(parameterName, description, ...
  parameterGroup, lockParam, dataType, numberOfDataDimensions, ...
  dataDimensions, data); %[]); %

nIndex = c3dServer.GetParameterIndex('POINT', 'FRAMES');
c3dServer.SetParameterValue(nIndex, 0, -1); %[]); %


%%
parameterGroup = 'TRIAL';
c3dServer.AddGroup(0, parameterGroup, '', '0')

parameterName = 'ACTUAL_END_FIELD';
description = '';
lockParam = '0';
dataType = int16(4); %2); %4; %-1; %-1 : Character, 1 : Byte, 2 : Integer, 4 : Floating point
numberOfDataDimensions = int16(2); %int16(1); %2; %2; %needs to be 2 for a string, since first dimension holds the length of each string
dataDimensions = int16([1 2]); %1; %{uint8(1)}; %[1 1]; %[1, numberOfMarkers]; % [lengthOfEachString, numberOfStrings]
data = int16([mod(numberOfFrames, (2^16 - 1)), floor(numberOfFrames / (2^16 - 1))]);

c3dServer.AddParameter(parameterName, description, ...
  parameterGroup, lockParam, dataType, numberOfDataDimensions, ...
  dataDimensions, data); %[]); %

parameterName = 'TEMP';
description = '';
lockParam = '0';
dataType = int16(4); %2); %4; %-1; %-1 : Character, 1 : Byte, 2 : Integer, 4 : Floating point
numberOfDataDimensions = int16(2); %int16(1); %2; %2; %needs to be 2 for a string, since first dimension holds the length of each string
dataDimensions = int16([1 2]); %1; %{uint8(1)}; %[1 1]; %[1, numberOfMarkers]; % [lengthOfEachString, numberOfStrings]
data = int16([0, 1]);

% c3dServer.AddParameter(parameterName, description, ...
%   parameterGroup, lockParam, dataType, numberOfDataDimensions, ...
%   dataDimensions, data); %[]); %


% nIndex = c3dServer.GetParameterIndex('TRIAL', 'ACTUAL_START_FIELD');
% nIndex
% % c3dServer.SetParameterName(nIndex, 'ACTUAL_START_FIELD');
% c3dServer.SetParameterValue(nIndex, 0, 0); %[]); %
% c3dServer.SetParameterValue(nIndex, 1, 1); %[]); %



%%


savec3d(c3dServer); %, 'output.c3d'); %
closec3d(c3dServer);

% phasespaceAndTreadmillStructure.markers = customPhasespaceLog;

phasespaceRecord = allMarkerData;
outputFilename = createMatOutputFilename(outputFilename);
save(outputFilename, 'phasespaceRecord');


end

