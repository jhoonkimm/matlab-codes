function [ phasespaceAndAnalogStructure ] = savePhasespaceAndAnalogDataInC3D( ...
  customPhasespaceLogName, tdmsFilename, tdmsTabNumber, c3dFilenameToLoad, outputFilename, varargin)
%SAVEPHASESPACEANDANALOGDATAINC3D requires a 32 bit matlab and an
%installation of c3dserver

addpath('./C3D');
addpath('../importTDMS');

% customPhasespaceLogName = ...
%   'C:\jrebula\myProjects\IMUGaitAnalysis\rawData\20110512-customCodeSyncTest\phasespaceCustomLog1-testTreadmillSync.txt';
% tdmsFilename = ...
%   'C:\jrebula\myProjects\LXL_5_10_11.tdms';
% tdmsTabNumber = 2; %1; %

verbose = 10;
shouldFillGapsInMarkerData = 1; %0; %
errorOnFindingNaNs = 0;
parameterDataToWrite = [];
forceReprocess = 0;
% c3dFilenameToLoad = [];
% outputFilename = 'output.c3d';
% if(0)
c3dUnitsAreInMeters = 0;
for i = 1:2:length(varargin)
  if (strcmp(varargin{i}, 'verbose'))
    verbose = varargin{i+1};
  end
  if (strcmp(varargin{i}, 'forceReprocess'))
    forceReprocess = varargin{i+1};
  end
  if (strcmp(varargin{i}, 'shouldFillGapsInMarkerData'))
    shouldFillGapsInMarkerData = varargin{i+1};
  end
  if (strcmp(varargin{i}, 'errorOnFindingNaNs'))
    errorOnFindingNaNs = varargin{i+1};
  end
  if (strcmp(varargin{i}, 'c3dUnitsAreInMeters'))
    c3dUnitsAreInMeters = varargin{i+1};
  end
  if (strcmp(varargin{i}, 'parameterDataToWrite'))
    parameterDataToWrite = varargin{i+1};
  end
  %   if (strcmp(varargin{i}, 'c3dFilenameToLoad'))
  %     c3dFilenameToLoad = varargin{i+1};
  %   end
  %   if (strcmp(varargin{i}, 'outputFilename'))
  %     outputFilename = varargin{i+1};
  %   end
end
% end

%%
% projectsDir = 'C:\Users\jrebula\myProjects\';

%%
% c3dFilenameToLoad = ...
%   ...[projectsDir 'CESR Project\Amputee Testing\Subject C2E_UM01_2011-05-25\PRES\C2E_UM01-2011-05-25-PRES-01.c3d'];
%   ...[projectsDir 'CESR Project\Amputee Testing\Subject C2E_UM01_2011-05-25\PRES\C2E_UM01-2011-05-25-fxn-joint-01.c3d'];
%   [projectsDir 'CESR Project\Amputee Testing\Subject C2E_UM01_2011-05-25\PRES\C2E_UM01-2011-05-25-PRES-02.c3d'];
% customPhasespaceLogName = ...
%   [projectsDir 'CESR Project\Amputee Testing\Subject C2E_UM01_2011-05-25\PRES\phasespaceCustomLog2.txt'];
%   ...[projectsDir 'CESR Project\Amputee Testing\Subject C2E_UM01_2011-05-25\PRES\phasespaceCustomLog1.txt'];
% tdmsFilename = ...
%   [projectsDir 'CESR Project\Amputee Testing\Subject C2E_UM01_2011-05-25\CESR Subject C2E_UM01_2011-05-25 PRES.tdms'];
% tdmsTabNumber = 2; %1; %


% dataDir = [projectsDir 'IMUGaitAnalysis\rawData\Rapture-InverseDynamicsTest\'];
%%
% outputFilename = 'rapture_1_25Combined.c3d';
% c3dFilenameToLoad = ...
%   [dataDir 'Tracked\take2-newexport.c3d'];
% customPhasespaceLogName = ...
%   [dataDir 'CPL\phasespaceCustomLog1.txt'];
% tdmsFilename = ...
%   [dataDir ' TDMS\rapture1.tdms'];
% tdmsTabNumber = 1; %2; %

% originalTrackedWalkTrial = loadPhasespaceRecord(c3dFilenameToLoad);
% for markerName = originalTrackedWalkTrial.markerNames
%   markerName{:}
%   sum(sum(isnan(originalTrackedWalkTrial.(markerName{:}))))
% end

%%
% outputFilename = 'functionalJointsCombined.c3d';
% customPhasespaceLogName = ...
%   [dataDir 'CPL\phasespaceCustomLog3.txt'];
% c3dFilenameToLoad = ...
%   [dataDir 'Tracked\take4-newexport.c3d']; %test-04.c3d']; %Tracked\take4.c3d']; %
% tdmsFilename = ...
%   [dataDir ' TDMS\functional joints.tdms'];
% tdmsTabNumber = 1; %2; %

% originalTrackedFunctionalTrial = loadPhasespaceRecord(c3dFilenameToLoad, 'forceReprocess', 1);
% for markerName = originalTrackedFunctionalTrial.markerNames, markerName{:}, sum(sum(isnan(originalTrackedFunctionalTrial.(markerName{:})))), end

%% 0.7m/s
% outputFilename = 'rapture_7Combined.c3d';
% c3dFilenameToLoad = ...
%   [dataDir 'Recap\test_05_trial_04.c3d'];
% customPhasespaceLogName = ...
%   [dataDir 'CPL\phasespaceCustomLog4.txt'];
% tdmsFilename = ...
%   [dataDir ' TDMS\rapture_7.tdms'];
% tdmsTabNumber = 1; %2; %

%% 1.6m/s
% outputFilename = 'rapture_1_8Combined.c3d';
% c3dFilenameToLoad = ...
%   [dataDir 'Recap\test_09_trial_08.c3d'];
% customPhasespaceLogName = ...
%   [dataDir 'CPL\phasespaceCustomLog8.txt'];
% tdmsFilename = ...
%   [dataDir ' TDMS\rapture1_6.tdms'];
% tdmsTabNumber = 1; %2; %

%%
if (verbose > 0)
  fprintf('loading in custom phasespace log...\n');
end
[customPhasespaceLog] = processCustomPhasespaceLog(customPhasespaceLogName, 'forceReprocess', forceReprocess);
markerNames = fields(customPhasespaceLog);
numberOfMarkers = length(markerNames);


%%
tdmsFile = importTDMSV3(tdmsFilename);

%%
groups = tdmsFile.group;
fprintf('tdms file has %g tabs\n', length(groups));
% for tdmsTabNumber = 1:length(groups)
lastGroup = groups(tdmsTabNumber);
channels = lastGroup.channel;
analogChannelNames = {};
for i = 1:length(channels)
  analogChannelNames{i} = strrep(channels(i).name, '/', '_');
  analogChannelNames{i} = ['A' analogChannelNames{i}((end-4):end)];
end
analogChannelData = {};
for i = 1:length(channels)
  analogChannelData{i} = channels(i).data;
end

%%
% parsePhasespaceTimeCode(timeCodeSignal, analogToMarkerRatio)

timeCodeSignal = analogChannelData{13}; %analogChannelData{14}; %
analogToMarkerRatio = 2; % since we are sampling at 960 = 480*2 Hz
inds = find(~(timeCodeSignal > 1));
ind = find(diff(inds) > 10 * 12 * analogToMarkerRatio, 1);
ia = inds(ind + 1) + 6 * analogToMarkerRatio; % double check the +1
indecesForTimingSignalData = ia + 12 * analogToMarkerRatio * [21:28, 11:18, 1:8]; %[1:8, 11:18, 21:28]; %
frameNumberBinary = timeCodeSignal(indecesForTimingSignalData); % double check this
% frameNumber = bin2dec(fliplr(regexprep(num2str(frameNumberBinary > 1), ' ', '')));
frameNumber = bin2dec(regexprep(num2str(frameNumberBinary > 1), ' ', ''));

frameNumberBinaryNext = timeCodeSignal(indecesForTimingSignalData + 480 * analogToMarkerRatio); % double check this
% frameNumberNext = bin2dec(fliplr(regexprep(num2str(frameNumberBinaryNext > 1), ' ', '')));
frameNumberNext = bin2dec(regexprep(num2str(frameNumberBinaryNext > 1), ' ', ''));

if (frameNumberNext - frameNumber ~= 480)
  fprintf('it doesn''t look like the phasespace frame rate is 480Hz, as seen from the timecode read in from the ADC (the tdms file)!');
end

% plot(timeCodeSignal)
% hold on;
% plot(indecesForTimingSignalData, timeCodeSignal(indecesForTimingSignalData), '*k')

firstTDMSTickInPhasespaceFrames = frameNumber - floor(ia / analogToMarkerRatio);
firstPhasespaceCustomFrame = customPhasespaceLog.frameNumber(1);

if (abs(firstPhasespaceCustomFrame - firstTDMSTickInPhasespaceFrames) > 100)
  %   fprintf('warning, there appears to be a difference of %g phasespace frames between the start of the tdms file and the custom phasespace log!\n', ...
  %     firstPhasespaceCustomFrame - firstTDMSTickInPhasespaceFrames);
  preposition = 'before';
  if (firstPhasespaceCustomFrame - firstTDMSTickInPhasespaceFrames < 0)
    preposition = 'after';
  end
  frameDiff = abs(firstPhasespaceCustomFrame - firstTDMSTickInPhasespaceFrames);
  fprintf('warning, the tdms file appears to start %g phasespace frames (%g s) %s the custom phasespace log!\n', ...
    frameDiff, frameDiff / 480, preposition);
end
% end
firstCommonFrame = max(firstTDMSTickInPhasespaceFrames, firstPhasespaceCustomFrame);

customPhasespaceLogSynced = customPhasespaceLog;
for i = 1:length(channels)
  analogChannelData{i} = ...
    channels(i).data(1 + analogToMarkerRatio * (firstCommonFrame - firstTDMSTickInPhasespaceFrames):end);
end
for markerNumber = 0:(numberOfMarkers - 1)
  markerName = markerNames{markerNumber + 1};
  if (strcmp(markerName, 'frameNumber'))
    continue;
  end
  markerData = customPhasespaceLogSynced.(markerName);
  customPhasespaceLogSynced.(markerName) = ...
    markerData(1 + (firstCommonFrame - firstPhasespaceCustomFrame):end, :);
end

%% trimming out first 6 samples of analog data, due to the stupid delay
% between the force measurements and the adc...
for i = 1:length(channels)
  analogChannelData{i} = analogChannelData{i}(6:end);
end

shortestLength = min([length(customPhasespaceLogSynced.(markerName)) * analogToMarkerRatio, ...
  length(analogChannelData{1})]);
shortestLength = shortestLength - mod(shortestLength, analogToMarkerRatio);

for i = 1:length(channels)
  analogChannelData{i} = ...
    analogChannelData{i}(1:shortestLength);
end
for markerNumber = 0:(numberOfMarkers - 1)
  markerName = markerNames{markerNumber + 1};
  if (strcmp(markerName, 'frameNumber'))
    continue;
  end
  markerData = customPhasespaceLogSynced.(markerName);
  customPhasespaceLogSynced.(markerName) = ...
    markerData(1 : shortestLength / analogToMarkerRatio, :);
end

if (diff([length(customPhasespaceLogSynced.(markerName)) * analogToMarkerRatio, ...
    length(analogChannelData{1})]))
  error('trimmed analog data length doesn''t match trimmed marker data length');
end

customPhasespaceLog = customPhasespaceLogSynced;

%%
% c3dServer = actxserver('C3DServer.C3D');

%%
% defaultC3dToLoad = 'C2E_A01 TM 120 A 02.c3d';



if (~isempty(c3dFilenameToLoad))
  trackedC3dTrial = loadPhasespaceRecord(c3dFilenameToLoad, 'forceReprocess', 1, ...
    'c3dUnitsAreInMeters', c3dUnitsAreInMeters, 'forceReprocess', forceReprocess);
  for markerNumber = 1:length(trackedC3dTrial.markerNames)
    markerName = trackedC3dTrial.markerNames{markerNumber};
    if(shouldFillGapsInMarkerData)
      if (sum(sum(isnan(trackedC3dTrial.(markerName)))) > 0)
        fprintf('warning, marker %s isn''t fully tracked in the c3d %s\n', markerName, c3dFilenameToLoad);
        
        [trackedData] = fillGapsInMarkerData(markerName, trackedC3dTrial.(markerName));
        trackedC3dTrial.(markerName) = trackedData;
      end
    end
  end
  
  %%
  %   currentMarkerNumber = 1;
  for markerNumber = 1:length(markerNames)
    markerName = markerNames{markerNumber};
    if (strcmp(markerName, 'frameNumber'))
      continue;
    end
    markerNameC3DCorrespondingToPhasespace{markerNumber} = sprintf('M%03g', str2double(markerName(2:end)));
    markerDataC3D = trackedC3dTrial.(markerNameC3DCorrespondingToPhasespace{markerNumber});
    %     markerDataC3D = get3dtarget(c3dServer, markerNameC3d, 0);
    %     allMarkersC3D.(markerNameC3d) = markerDataC3D;
    %     currentMarkerNumber = currentMarkerNumber + 1;
    %     break;
  end
  
  for markerNumber = 1:length(markerNames)
    markerName = markerNames{markerNumber};
    if (strcmp(markerName, 'frameNumber'))
      continue;
    end
    markerDataCustomLog = customPhasespaceLog.(markerName);
    markerDataC3D = trackedC3dTrial.(markerNameC3DCorrespondingToPhasespace{markerNumber});
    
    if(sum(isnan(markerDataCustomLog(:,1))) < length(markerDataCustomLog)*0.5)
      break;
    end
  end
  
  markerDataCustomLog = customPhasespaceLog.(markerName);
  markerDataC3D = [-markerDataC3D(:, 1), markerDataC3D(:, 3), markerDataC3D(:, 2)];
  
  markerDataC3D = markerDataC3D(1:min(length(markerDataC3D), length(markerDataCustomLog)), :);
  markerDataCustomLog = markerDataCustomLog(1:min(length(markerDataC3D), length(markerDataCustomLog)), :);
  
  %%
  %   [vals, lags] = xcorr(markerDataC3D(:, 1), markerDataCustomLog(:, 1));
  %   for markerIndexToUse = 1:size(markerDataCustomLog, 1)
  customLogDataToSyncWithTrackedC3D = markerDataCustomLog(:, 1);
  trackedDataToSyncWithCustomLog = markerDataC3D(:, 1);
  %     if((sum(isnan(customLogDataToSyncWithTrackedC3D))) < length(customLogDataToSyncWithTrackedC3D)/2)
  %       break;
  %     end
  %   end
  
  if (sum(isnan(customLogDataToSyncWithTrackedC3D)))
    isGoodVector = ~isnan(customLogDataToSyncWithTrackedC3D);
    indecesOfBeginningsOfGoodStretches = find(diff(isGoodVector) > 0.5) + 1;
    indecesOfEndsOfGoodStretches = find(diff(isGoodVector) < -0.5);
    %     lengthsOfGoodStretches = diff(indecesOfBeginningsOfGoodStretches);
    
    if (indecesOfEndsOfGoodStretches(1) < indecesOfBeginningsOfGoodStretches(1))
      indecesOfBeginningsOfGoodStretches = [1; indecesOfBeginningsOfGoodStretches];
    else
      % indecesOfEndsOfGoodStretches = [1; indecesOfEndsOfGoodStretches];
    end
    
    if (indecesOfEndsOfGoodStretches(end) < indecesOfBeginningsOfGoodStretches(end))
      indecesOfEndsOfGoodStretches = ...
        [indecesOfEndsOfGoodStretches; length(customLogDataToSyncWithTrackedC3D)];
    else
      %       indecesOfBeginningsOfGoodStretches = ...
      %         [indecesOfBeginningsOfGoodStretches; length(customLogDataToSyncWithTrackedC3D)];
    end
    
    lengthsOfGoodStretches = indecesOfEndsOfGoodStretches - indecesOfBeginningsOfGoodStretches;
    [lengthOfLongestStretch, indexOfMaxStretch] = max(lengthsOfGoodStretches);
    indecesOfGoodData = (5:(lengthOfLongestStretch - 5)) + indecesOfBeginningsOfGoodStretches(indexOfMaxStretch);
    if (sum(isnan(customLogDataToSyncWithTrackedC3D(indecesOfGoodData))))
      %       if (errorOnFindingNaNs)
      error('failed to find a period of custom log data without nans');
      %       end
    end
    %     cla;
    %     plot(customLogDataToSyncWithTrackedC3D(indecesOfGoodData))
    customLogDataToSyncWithTrackedC3D = customLogDataToSyncWithTrackedC3D(indecesOfGoodData);
    trackedDataToSyncWithCustomLog = trackedDataToSyncWithCustomLog(indecesOfGoodData);
  end
  
  %   firstIndexToUseInSync = floor(length(trackedDataToSyncWithCustomLog)/3);
  %   indecesToUseForSync = firstIndexToUseInSync:(firstIndexToUseInSync*2);
  firstIndexToUseInSync = 1;
  indecesToUseForSync = firstIndexToUseInSync:length(trackedDataToSyncWithCustomLog);
  
  %   %   figure;
  %   cla;
  %   plot(customLogDataToSyncWithTrackedC3D / 1000.0, 'LineWidth', 4);
  %   %   plot(diff(sum(markerDataCustomLog'.^2)), '--')
  %   hold on;
  %   %   plot(markerDataC3D, '--', 'LineWidth', 2);
  %   plot(indecesToUseForSync, trackedDataToSyncWithCustomLog(indecesToUseForSync), '--', 'LineWidth', 2, 'Color', 'k')
  %   %   plot(diff(sum(markerDataC3D'.^2)));
  %   %       legend({'markerDataC3D', 'markerDataCustomLog'});
  
  %%
  %   [vals, lags] = xcorr(markerDataC3D(:, 1), markerDataCustomLog(:, 1));
  
  [vals, lags] = xcorr( ...
    (customLogDataToSyncWithTrackedC3D - mean(customLogDataToSyncWithTrackedC3D))/1000, ...
    (trackedDataToSyncWithCustomLog(indecesToUseForSync) - mean(trackedDataToSyncWithCustomLog(indecesToUseForSync))), ...
    'none');
  
  %   [vals, lags] = xcorr( ...
  %     (customLogDataToSyncWithTrackedC3D)/1000, ...
  %     (trackedDataToSyncWithCustomLog(indecesToUseForSync)), ...
  %     'none');
  %   cla;
  %   plot(lags - firstIndexToUseInSync + 1, vals)
  [maxCorr, maxCorrInd] = max(vals);
  maxLag = lags(maxCorrInd);
  %   plot((1:length(markerDataC3D)) + maxLag, markerDataCustomLog(:, 1))
  
  range = max(max(markerDataC3D)) - min(min(markerDataC3D));
  if (range < 30)
    error(['it looks like the marker data from the c3d is either barely moving, or it''s not in mm, it''s in meters. Check this... ' ...
      'If the c3d is actually in meters, rerun this code adding to the function call: savePhasespaceAndAnalogDataInC3D(..., ''c3dUnitsAreInMeters'', 1)']);
  end
  
  differencesBetweenC3dMarkersAndCustomLog = markerDataCustomLog - markerDataC3D;
  if (max(max(abs(differencesBetweenC3dMarkersAndCustomLog))) > 1e-2)
    %     error('custom log and c3d log give different results, still have to implement time lag rectification, calculated max lag is %g\n', maxLag);
    fprintf('custom log and c3d log give different results, performing time lag rectification, calculated max lag is %g\n', maxLag);
    if (abs(maxLag) > 1000)
      error('the lag between the c3d marker data na dhte markerDataCustomLog is large, and the resulting curves don''t match well. Are you sure you are specifying c3d and customLogs that correspond to each other?');
    end
  end
  
  customLogIndecesToUse = 1:length(markerDataCustomLog);
  customLogIndecesToUse = customLogIndecesToUse - maxLag;
  
  customLogIndecesToUse(customLogIndecesToUse > length(markerDataCustomLog)) = [];
  customLogIndecesToUse(customLogIndecesToUse < 1) = [];
  
  %% swap in the c3d data for the custom log data:
  for markerNumber = 1:numberOfMarkers
    markerName = markerNames{markerNumber};
    if (strcmp(markerName, 'frameNumber'))
      continue;
    end
    markerNameC3DCorrespondingToPhasespace{markerNumber} = sprintf('M%03g', str2double(markerName(2:end)));
    markerDataC3D = trackedC3dTrial.(markerNameC3DCorrespondingToPhasespace{markerNumber});
    markerDataC3D = markerDataC3D(customLogIndecesToUse, :);
    
    %     markerNameC3d = sprintf('M%03g', str2double(markerName(2:end)));
    %     markerNameC3d
    %     try
    %       markerDataC3D = get3dtarget(c3dServer, markerNameC3d, 0);
    %     markerDataC3D = allMarkersC3D.(markerNameC3d);
    %     catch e
    %       e
    %     end
    markerDataC3D = [-markerDataC3D(:, 1), markerDataC3D(:, 3), markerDataC3D(:, 2)];
    
    markerDataCustomLog = customPhasespaceLog.(markerName);
    
    markerDataC3D = markerDataC3D(1:min(length(markerDataC3D), length(markerDataCustomLog)), :);
    markerDataCustomLog = markerDataCustomLog(1:min(length(markerDataC3D), length(markerDataCustomLog)), :);
    
    if (sum(sum(isnan(markerDataC3D))))
      if (errorOnFindingNaNs)
        error('shouldn''t replace marker data with nans!');
      end
    end
    
    customPhasespaceLog.(markerName) = markerDataC3D * 1000.0;
    %     cla;
    %     plot(markerDataC3D * 1000.0)
    %     hold on
    %     plot(customPhasespaceLog.(markerName), '--k');
    
    %     fprintf('overwriting untracked marker %s with tracked\n', markerName);
    %     currentMarkerNumber = currentMarkerNumber + 1;
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
videoFrameRate = 480;
numberOfFrames = length(customPhasespaceLog.(markerNames{1}));
analogToVideoRatio = analogToMarkerRatio;
numberOfAnalogChannels = length(channels);
filetype = 3; %1; %
dataType = 2; %
pointDataScaleFactor = 1.0;

createc3d(c3dServer, name, videoFrameRate, numberOfFrames, numberOfMarkers - 1, ...
  analogToVideoRatio, numberOfAnalogChannels, filetype, dataType, pointDataScaleFactor);

%%
% first rip out a bunch of data from the template c3d:

if (0)
  analogChannelsTemplate = c3dServer.getAnalogChannels();
  numberOfAnalogChannelsTemplate = length(fields(analogChannelsTemplate));
  for analogChannelIndex = (numberOfAnalogChannelsTemplate - 1):0
    try
      c3dServer.DeleteAnalogChannelWithParam(0);
      %   c3dServer.DeleteAnalogChannel(0);
    catch e
      e
      break;
    end
  end
  for analogChannelIndex = (numberOfAnalogChannelsTemplate - 1):0
    c3dServer.DeleteAnalogChannelWithParam(0);
  end
  numberOfMarkersTemplate = c3dServer.GetNumber3DPoints;
  for analogChannelIndex = 0:(numberOfMarkersTemplate - 1)
    c3dServer.DeleteMarker(0);
  end
  
  firstFrame = c3dServer.GetVideoFrame(0);
  numFrames = c3dServer.GetVideoFrame(1) - firstFrame + 1;
  c3dServer.DeleteFrames(firstFrame, numFrames);
  
  % c3dServer.SetAVRatio(2)
  
  savec3d(c3dServer, 'template.c3d');
  return
end

%% write analogs:


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
      dataDimensions, data); %[]); %
    %% Now Save the Modified file with the new Offset regions:
    %   c3dServer.AddParameter(parameterName,'Array of Baseline Indices for Multiple Force Platforms', 'FORCE_PLATFORM', '0', ...
    %     int16(2), int16(2), int16([2 2]), int16([zeroindsfp1{qq} zeroindsfp2{qq}]));
  end
end


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
nIndex = c3dServer.GetParameterIndex('ANALOG', 'LABELS');
for analogChannelNumber = 0:(length(analogChannelNames) - 1)
  analogChannelName = analogChannelNames{analogChannelNumber + 1}; %(((end-5):end)); %
  %   fprintf('setting analog channel name: %s\n', analogChannelName);
  c3dServer.SetParameterValue(nIndex, analogChannelNumber, analogChannelName);
end

%% write marker data

markerIndex = 0;
for markerNumber = 1:numberOfMarkers
  markerName = markerNames{markerNumber};
  if (strcmp(markerName, 'frameNumber'))
    continue;
  end
  %   markerIndex
  %   markerName
  markerData = customPhasespaceLog.(markerName);
  
  if (sum(sum(isnan(markerData))))
    if (errorOnFindingNaNs)
      error('shoulnd''t output nan data from filled in c3d!');
    end
  end
  %   fprintf('correctly removed nans from marker %s\n', markerName);
  
  
%   if (shouldFillGapsInMarkerData)
%     markerData = fillGapsInMarkerData(markerName, markerData);
%   end
  
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
      % we want to do a transformation:
      % y->z
      % -z->y
      % x->x
      if (channelToWrite == 0)
        channelToWrite = 0;
      elseif (channelToWrite == 1)
        channelToWrite = 2;
      elseif (channelToWrite == 2)
        channelToWrite = 1;
        dataToWrite = -1.0 * dataToWrite;
      end
      
      %       if(fillMarkerGapsWithSplines)
      %         goodIndeces = ~isnan(dataToWrite);
      %         allIndeces = 1:length(dataToWrite);
      %         dataToWrite = spline(allIndeces(goodIndeces), dataToWrite(goodIndeces), allIndeces);
      %       end
      
      %       sum(isnan(dataToWrite))
      dataToWrite(isnan(dataToWrite)) = 0;
      if(c3dUnitsAreInMeters)
        markerDataOutputScale = 0.001;
      else
        markerDataOutputScale = 1;
      end
      
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

savec3d(c3dServer); %, 'output.c3d'); %
closec3d(c3dServer);

phasespaceAndAnalogStructure.markers = customPhasespaceLog;

for i = 1:length(analogChannelNames)
  analogChannelName = analogChannelNames{i}; %(((end-5):end)); %
  phasespaceAndAnalogStructure.analogChannels.(analogChannelName) = ...
    analogChannelData{i};
end

%%
a = loadPhasespaceRecord(outputFilename, 'forceReprocess', 1);
totalBad = 0;
for markerName = a.markerNames
  %   markerName{:}
  totalBad = totalBad + sum(sum(isnan(a.(markerName{:}))));
end
if (totalBad > 0)
  if (errorOnFindingNaNs)
    error('reloaded processed c3d, and found missing markers, not fully tracked! This is a programming error, check the code!');
  end
end


end


