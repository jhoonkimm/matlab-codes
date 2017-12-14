function mergeC3D(varargin)


errorOnFindingNaNs = 0;
phasespaceSampleFrequency=480;
c3dUnitsAreInMeters = 1;
forceReprocess = 0;
parameterDataToWrite = [];

numFiles=length(varargin)-1;
inputFiles=varargin(1:end-1);
outputFilename=varargin{end};


%% load first file, build initial struct
c3dFilenameToLoad = varargin{1};
trackedC3dTrial = loadPhasespaceRecord(c3dFilenameToLoad, ...
        'c3dUnitsAreInMeters', c3dUnitsAreInMeters, 'forceReprocess', forceReprocess);
    
if (~isfield(trackedC3dTrial, 'markerNames'))
    trackedC3dTrial.markerNames = fieldnames(trackedC3dTrial);
end

rootMarkerNames = trackedC3dTrial.markerNames(:);
rootNumberOfMarkers = length(rootMarkerNames);
rootTrackedC3D = trackedC3dTrial;
numFiles = numFiles - 1;
varargin(1) = [];
%% load following files, append to struct

while numFiles > 0
  c3dFilenameToLoad = varargin{1};
  trackedC3dTrial = loadPhasespaceRecord(c3dFilenameToLoad, ...
          'c3dUnitsAreInMeters', c3dUnitsAreInMeters, 'forceReprocess', forceReprocess);

  if (~isfield(trackedC3dTrial, 'markerNames'))
      trackedC3dTrial.markerNames = fieldnames(trackedC3dTrial);
  end
  
  if isequal(rootMarkerNames, trackedC3dTrial.markerNames(:)) %markers match
    for m=1:rootNumberOfMarkers
      rootTrackedC3D.(rootMarkerNames{m}) = vertcat(rootTrackedC3D.(rootMarkerNames{m}),trackedC3dTrial.(rootMarkerNames{m}));
    end
  end
  
  numFiles = numFiles - 1;
  varargin(1) = [];
end

numberOfMarkers = rootNumberOfMarkers;
markerNames = rootMarkerNames;
%% Write final file
c3dServer = actxserver('C3DServer.C3D');
name = outputFilename; %'output.c3d';
videoFrameRate = phasespaceSampleFrequency;
numberOfFrames = length(rootTrackedC3D.(markerNames{find(~strcmp(markerNames, 'frameNumber'),1)}));
analogToVideoRatio = 1;
numberOfAnalogChannels = 0;
filetype = 3; %1; %
dataType = 2; %
pointDataScaleFactor = 1.0;

createc3d(c3dServer, name, videoFrameRate, numberOfFrames, numberOfMarkers, ... - 1, ...
    analogToVideoRatio, numberOfAnalogChannels, filetype, dataType, pointDataScaleFactor);



%% write analogs:

% nIndex = c3dServer.GetParameterIndex('ANALOG', 'OFFSET');
% c3dServer.SetParameterType(nIndex, 4);
% analogIndex = 0;
% for analogNumber = 1:numberOfAnalogChannels
%     c3dServer.SetParameterValue(nIndex, analogIndex, single(0.0));
%     analogIndex = analogIndex + 1;
% end
% 
% 
% for channelNumber = 0:(numberOfAnalogChannels - 1)
%     %   fprintf('writing channel %g\n', channelNumber);
%     
%     dataToWrite = zeros(1, numberOfFrames * 2);
%     dataToWrite(1:length(single(analogChannelData{channelNumber + 1}))) = ...
%         single(analogChannelData{channelNumber + 1});
%     
%     %   subFrameNums = 1:(numberOfFrames * analogToVideoRatio);
%     c3dServer.SetAnalogDataEx(channelNumber, 1, single(dataToWrite + 0)); %
%     % sin(subFrameNums/10)));
%     
%     %   frameNums = 1:(numberOfFrames); for frameNum = frameNums
%     %     for subFrameNum = 1:analogToVideoRatio
%     %       c3dServer.SetAnalogData(channelNumber, frameNum, subFrameNum,
%     %       sin(frameNum/10));
%     %     end
%     %   end
% end

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
%% set the length of the trial explicitly, so that we can have trials longer than 2^16 frames
parameterGroup = 'TRIAL'; c3dServer.AddGroup (0, parameterGroup, '', '0');

parameterName = 'ACTUAL_START_FIELD'; data = int16([1, 0]);
description = '';
lockParam = '0';
dataType = int16(2); %2); %4; %-1; %-1 : Character, 1 : Byte, 2 :
% Integer, 4 : Floating point
numberOfDataDimensions = int16(2); %int16(2);
% %2; %2; %needs to be 2 for a string, since first dimension holds the
% length of each string
dataDimensions = int16([1 2]); %1; %{uint8(1)}; %[1
% 1]; %[1, numberOfMarkers]; % [lengthOfEachString, numberOfStrings]
%
c3dServer.AddParameter(parameterName, description, ...
  parameterGroup, lockParam, dataType, numberOfDataDimensions, ...
  dataDimensions, data); %[]); %
%
c3dServer.SetStrictParameterChecking(0);
parameterName = 'ACTUAL_END_FIELD';
data = typecast(uint32(numberOfFrames),'uint16');
%
c3dServer.AddParameter(parameterName, description, ...
  parameterGroup, lockParam, dataType, numberOfDataDimensions, ...
  dataDimensions, data); %[]); %
%
%
% % nIndex = c3dServer.GetParameterIndex('TRIAL', 'ACTUAL_START_FIELD'); %
% c3dServer.SetParameterType(nIndex, 2); % analogIndex = 0; % for
% analogNumber = 1:numberOfAnalogChannels %
% c3dServer.SetParameterValue(nIndex, analogIndex, single(0.0)); %
% analogIndex = analogIndex + 1; % end

% parameterGroup = 'POINT'; % c3dServer.AddGroup (0, parameterGroup, '',
% '0');
%
% parameterName = 'LONG_FRAMES'; description = ''; lockParam = '0';
% dataType = int16(2); %2); %4; %-1; %-1 : Character, 1 : Byte, 2 :
% Integer, 4 : Floating point numberOfDataDimensions = int16(2); %int16(2);
% %2; %2; %needs to be 2 for a string, since first dimension holds the
% length of each string dataDimensions = int16([1 1]); %1; %{uint8(1)}; %[1
% 1]; %[1, numberOfMarkers]; % [lengthOfEachString, numberOfStrings] data =
% numberOfFrames; %int16([mod(numberOfFrames, (2^16 - 1)),
% floor(numberOfFrames / (2^16 - 1))]);
%
% c3dServer.AddParameter(parameterName, description, ...
%   parameterGroup, lockParam, dataType, numberOfDataDimensions, ...
%   dataDimensions, data); %[]); %


%%
% parameterGroup = 'POINT'; % c3dServer.AddGroup (0, parameterGroup, '',
% '0');
%
% parameterName = 'LONG_FRAMES'; description = ''; lockParam = '0';
% dataType = int16(4); %2); %4; %-1; %-1 : Character, 1 : Byte, 2 :
% Integer, 4 : Floating point numberOfDataDimensions = int16(2); %int16(2);
% %2; %2; %needs to be 2 for a string, since first dimension holds the
% length of each string dataDimensions = int16([1 1]); %1; %{uint8(1)}; %[1
% 1]; %[1, numberOfMarkers]; % [lengthOfEachString, numberOfStrings] data =
% single([numberOfFrames 0]); %int16([mod(numberOfFrames, (2^16 - 1)),
% floor(numberOfFrames / (2^16 - 1))]);
%
% c3dServer.AddParameter(parameterName, description, ...
%   parameterGroup, lockParam, dataType, numberOfDataDimensions, ...
%   dataDimensions, data); %[]); %

for i = 1:length(parameterDataToWrite)
    parameterGroup = 'EXTRA';
    if (i == 1)
        c3dServer.AddGroup (0, parameterGroup, 'extra parameters added by matlab exporter', '0');
    end
    parameterName = parameterDataToWrite{i}{1}; %'ANGLE';
    description = '';
    lockParam = '0';
    dataType = int16(4); %2); %4; %-1; %-1 : Character, 1 : Byte, 2 : Integer, 4 : Floating point
    numberOfDataDimensions = int16(2); %int16(2); %2; %2; %needs to be 2 for a string, since first dimension holds the length of each string
    dataDimensions = int16([1 1]); %1; %{uint8(1)}; %[1 1]; %[1, numberOfMarkers]; % [lengthOfEachString, numberOfStrings]
    data = single([parameterDataToWrite{i}{2} 0]);
    
    addParametersToC3D = 1; %0; %
    if (addParametersToC3D)
        c3dServer.AddParameter(parameterName, description, ...
            parameterGroup, lockParam, dataType, numberOfDataDimensions, ...
            dataDimensions, data);
    end
end



%% name markers
nIndex = c3dServer.GetParameterIndex('POINT', 'LABELS');
% nItems = c3dServer.GetParameterLength(nIndex); target_name =
% c3dServer.GetParameterValue(nIndex, 0);
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
% nIndex = c3dServer.GetParameterIndex('ANALOG', 'LABELS');
% for analogChannelNumber = 0:(length(analogChannelNames) - 1)
%     analogChannelName = analogChannelNames{analogChannelNumber + 1}; %(((end-5):end)); %
%     %   fprintf('setting analog channel name: %s\n', analogChannelName);
%     c3dServer.SetParameterValue(nIndex, analogChannelNumber, analogChannelName);
% end

%% write marker data

markerIndex = 0;
for markerNumber = 1:numberOfMarkers
    markerName = markerNames{markerNumber};
    if (strcmp(markerName, 'frameNumber'))
        continue;
    end
    
    markerData = rootTrackedC3D.(markerName);
    
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
            % we want to do a transformation: y->z -z->y x->x 
            %% not for this app? %%
            if (channelToWrite == 0)
                channelToWrite = 0;
            elseif (channelToWrite == 1)
                channelToWrite = 1;
            elseif (channelToWrite == 2)
                channelToWrite = 2;
%                 dataToWrite = -1.0 * dataToWrite;
            end
            
            dataToWrite(isnan(dataToWrite)) = 0;
            
            markerDataOutputScale = 1;
            
            retVal = c3dServer.SetPointDataEx(markerIndex, channelToWrite, 1, single(dataToWrite * markerDataOutputScale));
        end
        
    end
    markerIndex = markerIndex + 1;
end


%%
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


savec3d(c3dServer); %, 'output.c3d'); %
closec3d(c3dServer);
