function [ phasespaceAndTreadmillStructure ] = savePhasespaceAndTreadmillDataInC3D2017(...
    customPhasespaceLogName, savedTreadmillDataFilename, c3dFilenameToLoad, outputFilename, varargin)
% savePhasespaceAndTreadmillDataInC3D2017 is a unified program to
% integrate Phasespace C3D, customlog, and treadmill data (either Lanxing's
% Matlab output or Labview TDMS output).
%
% Flags: 'verbose'
%   Default is true.
% 'forceReprocess'
%   Default is true. Use updated processCustomPhasespaceLog so that this
%   doesn't take forever.
% 'shouldFillGapsInMarkerData'
%   Default is true.
% 'errorOnFindingNaNs'
%   Default is false.
% 'c3dUnitsAreInMeters'
%   Default is false. Current standard output may require this flag.
% 'parameterDataToWrite' 'treadmillDataFromMatlab'
%   Default is true. Set to 0 to process TDMS files. Not actively
%   supported.
% 'tdmsTabNumber'
%   Default is 1. If using this, make sure the previous option is set to 0.
% 'phasespaceSampleFrequency'
%   Default is 480 Hz. Should be able to handle 240. Other numbers not
%   currently supported.
% 'ignoreStillTrials'
%   Default is false. Used to allow still trials to not error.
% 'truncateFrames'
%   Default is 0. Use to cut trials to length. Input a treadmill data (not
%   phasespace frame) file sample number to cut all data before, or two
%   numbers in an array to cut data before and after.
% 'c3dMarkerNamesHaveThreePaddedZeros',
%   Default is true, if false, this assumes that the tracked c3d files and
%   the custom phasespace log have the same marker nframes (or tracked is the
%   all uppercase version of the custom log), after all renaming has been
%   done.
% 'perturbationFile'
%   Default is empty. Path to perturbation file in xls format. TDMS
%   compatibility to be added later. 2015.07.22
%
% This program requires a 32 bit matlab and an installation of c3dserver
% (32 bit version). 64 bit beta version, however, has been known to work
%
% Last updated 2013.07.13 Minor updated 2014.08.05 by Amy Wu for marker
% renaming

% addpath('./C3D');

%% Default variables and conditions
verbose = 1;
shouldFillGapsInMarkerData = 1;
errorOnFindingNaNs = 0;
parameterDataToWrite = [];
forceReprocess = 1;
c3dUnitsAreInMeters = 0;
treadmillDataFromMatlab = 1;
tdmsTabNumber = 1;
phasespaceSampleFrequency = 480;
ignoreStillTrials = 0;
truncateFrames = 0;
markerRenaming = '';
c3dMarkerNamesHaveThreePaddedZeros = 0;
renameMarkersFromC3D = 0;
perturbationFile = [];
ForcePlateCorrection = 0;
PolyFitOrder = 2;
ForceFilter = 0;
ForceFilterFreq = 40; %Hz
mass = []; 
SwitchForcePlates = 0;
activity = 'Walk'; %Options are currently 'Run' and 'Walk.  This is only used if ForcePlateCorrection = 1

%% Optional input arguments
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
        case 'shouldFillGapsInMarkerData'
            shouldFillGapsInMarkerData = val;
        case 'errorOnFindingNaNs'
            errorOnFindingNaNs = val;
        case 'c3dUnitsAreInMeters'
            c3dUnitsAreInMeters = val;
        case 'parameterDataToWrite'
            parameterDataToWrite = val;
        case 'treadmillDataFromMatlab'
            treadmillDataFromMatlab = val;
        case 'tdmsTabNumber'
            tdmsTabNumber = val;
        case 'phasespaceSampleFrequency'
            phasespaceSampleFrequency = val;
            if mod(phasespaceSampleFrequency,240)
                error('phasespaceSampleFrequency invalid. Needs to be 240 or 480.')
            end
        case 'ignoreStillTrials';
            ignoreStillTrials = val;
        case 'truncateFrames'
            truncateFrames = val;
        case 'markerRenaming'
            markerRenaming = val;
        case 'renameMarkersFromC3D'
            renameMarkersFromC3D = val;
        case 'c3dMarkerNamesHaveThreePaddedZeros'
            c3dMarkerNamesHaveThreePaddedZeros = val;
        case 'perturbationFile'
            perturbationFile = val;
        case 'ForcePlateCorrection'
            ForcePlateCorrection = val;
        case 'PolyFitOrder'
            PolyFitOrder = val;
        case 'activity'
            activity = val;
        case 'ForceFilter'
            ForceFilter = val;
        case 'ForceFilterFreq'
            ForceFilterFreq = val;
        case 'SwitchForcePlates'
            SwitchForcePlates = val;
        case 'mass'
            mass = val;
        otherwise
            %       error('\nError: unrecognized parameter ''%s''\n',opt)
    end
end

%% Load custom log
if (verbose > 0) % barely used, should reconsider
    fprintf('loading in custom phasespace log...\n');
end
[customPhasespaceLog] = processCustomPhasespaceLog(customPhasespaceLogName, 'forceReprocess', forceReprocess,'markerRenaming',markerRenaming);

% rename markers in custom log:
markerNamesOriginal = fields(customPhasespaceLog);
for i = 1:length(markerNamesOriginal)
    renamed = 0;
    for j = 1:length(markerRenaming)
        if (strcmp(markerRenaming{j}{1}, markerNamesOriginal{i}))
            customPhasespaceLog2.([markerRenaming{j}{2}]) = customPhasespaceLog.(markerRenaming{j}{1});
            renamed = 1;
            break;
        end
    end
    if (~renamed)
        [customPhasespaceLog2.(markerNamesOriginal{i})] = customPhasespaceLog.(markerNamesOriginal{i});
    end
end
customPhasespaceLog = customPhasespaceLog2;

allFields = fields(customPhasespaceLog);
if isempty(markerRenaming)
    markerNames = allFields(strncmp(allFields,'M',1)); % John's original line
else
    markerNames = allFields(~strncmp(allFields,'frameNumber',3)); % Amy's mod for marker renaming
end
numberOfMarkers = length(markerNames);

if treadmillDataFromMatlab
    %% Load treadmill data from Matlab
    load(savedTreadmillDataFilename)
    
    %%
    analogChannelNames = {};
    for i = 1:size(samples, 2)
        analogChannelNames{i} = ['A' sprintf('%02g', (i-1))];
    end
    analogChannelData = {};
    for i = 1:size(samples, 2)
        analogChannelData{i} = samples(:, i)';
    end
    
    % reorder Matlab analog data to match TDMS, avoid needing second V3D
    % pipeline (fp1{x y z mx my mz} fp2{x y z mx my mz} timecode unusused)
    analogChannelData(7:14) = analogChannelData([9:14 8 7]);
else
    %% Load treadmill data from TDMS
    addpath('./importTDMS');
    tdmsFile = importTDMSV3(tdmsFilename);
    
    %%
    groups = tdmsFile.group;
    fprintf('tdms file has %g tabs\n', length(groups));
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
    
end

%% check timing
% % if treadmillDataFromMatlab
%   timeCodeSignal = analogChannelData{8};
% else
timeCodeSignal = analogChannelData{13};
% end

% read time code frame number (3 bytes chunks) from treadmill ADC, then at
% 1 second at treadmill ADC sampling rate later to verify phasespace
% counting frame rate. ** note that this is usually always 480, regardless
% of sampling frequency (which can be 240)
analogToMarkerRatio = 960/phasespaceSampleFrequency; % since we are sampling at 960 = phasespaceSampleFrequency*analogToMarkerRatio Hz
ratioFor480 = 960/480;
inds = find(~(timeCodeSignal > 1));
ind = find(diff(inds) > 10 * 12 * ratioFor480, 1); % find a start point after an appropriate pause (low signal) in time code
ia = inds(ind + 1) + 6 * ratioFor480; % skip a chunk (header?)
indecesForTimingSignalData = ia + 12 * ratioFor480 * [21:28, 11:18, 1:8]; % get 3 bytes of frame number in correct endian order
frameNumberBinary = timeCodeSignal(indecesForTimingSignalData);
frameNumber = bin2dec(regexprep(num2str(frameNumberBinary > 1), ' ', ''));

frameNumberBinaryNext = timeCodeSignal(indecesForTimingSignalData + 960); % double check this
frameNumberNext = bin2dec(regexprep(num2str(frameNumberBinaryNext > 1), ' ', ''));

if (frameNumberNext - frameNumber ~= 480)
    fprintf('it doesn''t look like the phasespace frame count rate is %d Hz, as seen from the timecode read in from the ADC (the tdms/mat file)! This program will probably not function properly.\n',phasespaceSampleFrequency);
end

% real phasespace frequency check
phasespaceFrameDiff = unique(diff(customPhasespaceLog.frameNumber));
expectedPhasespaceDiff = 480/phasespaceSampleFrequency;

if expectedPhasespaceDiff ~= phasespaceFrameDiff
    seenFrequency = 480/phasespaceFrameDiff;
    fprintf('Customlog seems to see a sample frequency of %d Hz instead of expected %d Hz\n', seenFrequency, phasespaceSampleFrequency)
    error('Check data or change phasespaceSampleFrequency input')
end

%% find and drop data before first common frame of treadmill and phasepace custom log
firstTreadmillTickInPhasespaceFrames = frameNumber - ia / ratioFor480;
firstPhasespaceCustomFrame = customPhasespaceLog.frameNumber(1);

if (abs(firstPhasespaceCustomFrame - firstTreadmillTickInPhasespaceFrames) > 100)
    preposition = 'before';
    if (firstPhasespaceCustomFrame - firstTreadmillTickInPhasespaceFrames < 0)
        preposition = 'after';
    end
    frameDiff = firstPhasespaceCustomFrame - firstTreadmillTickInPhasespaceFrames;
    fprintf('warning, the treadmill file appears to start %g phasespace frames (%g s) %s the custom phasespace log!\n', ...
        frameDiff, frameDiff / phasespaceSampleFrequency, preposition);
end
% end

firstCommonFrame = max(firstTreadmillTickInPhasespaceFrames, firstPhasespaceCustomFrame);

% At 240 Hz, phasespace frames increments by 2's, so odd frames don't
% exist. If first common frame is not a captured phasespace frame (can
% happen if the treadmill lags phasespace), find first existing phasespace
% frame after. This should also work for 120 Hz, theoretically.
firstCommonFrame = customPhasespaceLog.frameNumber(find(customPhasespaceLog.frameNumber>=firstCommonFrame,1));

analogFramesStart = 1 + ratioFor480 * (firstCommonFrame - firstTreadmillTickInPhasespaceFrames);
markerFramesStart = 1 + (firstCommonFrame - firstPhasespaceCustomFrame)/phasespaceFrameDiff;

customPhasespaceLogSynced = customPhasespaceLog;
for i = 1:length(analogChannelData)
    analogChannelData{i} = ...
        analogChannelData{i}(analogFramesStart:end);
end
for markerNumber = 1:numberOfMarkers
    markerName = markerNames{markerNumber};
    if (strcmp(markerName, 'frameNumber'))
        continue;
    end
    markerData = customPhasespaceLogSynced.(markerName);
    customPhasespaceLogSynced.(markerName) = ...
        markerData(markerFramesStart:end, :);
end

%% Obsolete timing correction due to Justin's correction code?
% trimming out first 6 samples of analog data, due to the stupid delay
% between the force measurements and the adc... for i =
% 1:length(analogChannelData)
% 	analogChannelData{i} = analogChannelData{i}(6:end);
% end

%% take shorter of customlog or treadmill data
shortestLength = min([length(customPhasespaceLogSynced.(markerNames{find(~strcmp(markerNames,'frameNumber'),1)})) * analogToMarkerRatio, ...
    length(analogChannelData{1})]);
shortestLength = shortestLength - mod(shortestLength, analogToMarkerRatio);

for i = 1:length(analogChannelData)
    analogChannelData{i} = ...
        analogChannelData{i}(1:shortestLength);
end
for markerNumber = 1:numberOfMarkers
    markerName = markerNames{markerNumber};
    if (strcmp(markerName, 'frameNumber'))
        continue;
    end
    markerData = customPhasespaceLogSynced.(markerName);
    customPhasespaceLogSynced.(markerName) = ...
        markerData(1 : shortestLength / analogToMarkerRatio, :);
end

if (diff([length(customPhasespaceLogSynced.(markerName)) * analogToMarkerRatio, ...
        length(analogChannelData{1})]))
    error('trimmed analog data length doesn''t match trimmed marker data length'); %shouldn't really happen
end

% truncate treadmill data if called
if any(truncateFrames) > 0
    analogTruncate = truncateFrames(1) - analogFramesStart;
    analogTruncate = ceil(analogTruncate/analogToMarkerRatio)*analogToMarkerRatio+1;
    if analogTruncate > length(analogChannelData{i})
        error('Truncate start begins after last available treadmill data point');
    elseif analogTruncate < 1
        disp('Warning: Truncate begins before treadmill and marker data sync begins. Using earliest available data point');
        analogTruncate = 1;
    end
    if length(truncateFrames)>1
        analogTruncate(2) = truncateFrames(2) - (analogFramesStart-1);
        analogTruncate(2) = analogTruncate(2) - mod(analogTruncate(2),analogToMarkerRatio);
        if analogTruncate(2) > length(analogChannelData{i})
            analogTruncate(2) = length(analogChannelData{i});
            disp('Warning: Truncate end is after the last available treadmill data point. Using last data point instead.')
        elseif analogTruncate(2) < 1
            error('Truncate end is before any good data starts.');
        end
    end
    for i = 1:length(analogChannelData)
        if length(analogTruncate)>1
            analogChannelData{i} = analogChannelData{i}(analogTruncate(1):analogTruncate(2));
        else
            analogChannelData{i} = analogChannelData{i}(analogTruncate(1):end);
        end
    end
end

customPhasespaceLog = customPhasespaceLogSynced;

%% open c3d file, optionally fill in gaps with splines
if (~isempty(c3dFilenameToLoad))
    trackedC3dTrial = loadPhasespaceRecord(c3dFilenameToLoad, ...
        'c3dUnitsAreInMeters', c3dUnitsAreInMeters, 'forceReprocess', forceReprocess);
    
    if (~isfield(trackedC3dTrial, 'markerNames'))
        trackedC3dTrial.markerNames = fieldnames(trackedC3dTrial);
    end
    
    % rename markers
    if (renameMarkersFromC3D)
        markerNamesOriginal = trackedC3dTrial.markerNames;
        trackedC3dTrial = rmfield(trackedC3dTrial, 'markerNames');
        
        %   trackedC3dTrial2 = trackedC3dTrial;
        for i = 1:length(markerNamesOriginal)
            renamed = 0;
            for j = 1:length(markerRenaming)
                if (strcmp(markerRenaming{j}{1}, markerNamesOriginal{i}))
                    [trackedC3dTrial2.(markerRenaming{j}{2})] = trackedC3dTrial.(markerRenaming{j}{1});
                    renamed = 1;
                    break;
                end
            end
            if (~renamed)
                [trackedC3dTrial2.(markerNamesOriginal{i})] = trackedC3dTrial.(markerNamesOriginal{i});
            end
        end
        trackedC3dTrial = trackedC3dTrial2;
        trackedC3dTrial.markerNames = fieldnames(trackedC3dTrial);
    end
    
    % use the c3d marker names now
    markerNames = trackedC3dTrial.markerNames(:);
    numberOfMarkers = length(markerNames);
    
    %   if(shouldFillGapsInMarkerData)
    %     for markerNumber = 1:length(trackedC3dTrial.markerNames)
    %       markerName = trackedC3dTrial.markerNames{markerNumber}; if
    %       (sum(sum(isnan(trackedC3dTrial.(markerName)))) > 0)
    %         fprintf('warning, marker %s isn''t fully tracked in the c3d
    %         %s\n', markerName, c3dFilenameToLoad);
    %
    %         [trackedData] = fillGapsInMarkerData(markerName,
    %         trackedC3dTrial.(markerName)); trackedC3dTrial.(markerName) =
    %         trackedData;
    %       end
    %     end
    %   end
    
    %% Pick a marker in the custom log that doesn't have a lot of nan's to check sync with
    for markerNumber = 1:length(markerNames)
        markerName = markerNames{markerNumber};
        if (strcmp(markerName, 'frameNumber'))
            continue;
        end
        try
            markerDataC3D = trackedC3dTrial.(markerName);
        catch e
            markerDataC3D = trackedC3dTrial.(upper(markerName));
        end
        try
            markerDataCustomLog = customPhasespaceLog.(markerName);
        catch e
            markerNames = fieldnames(customPhasespaceLog);
            markerName = markerNames(find(~cellfun(@isempty, strfind(upper(markerNames), markerName))));
            markerName = markerName{:};
            markerDataCustomLog = customPhasespaceLog.(markerName);
        end
        if(sum(isnan(markerDataCustomLog(:,1))) < length(markerDataCustomLog)*0.5)
            if sum(isnan(markerDataC3D(:,1))) == 0 % we need mostly good customlog and fully filled c3d, or else xcorr kills us later
                break;
            end
        end
    end
    
    try
        markerDataCustomLog = customPhasespaceLog.(markerName);
    catch e
        e
    end
    markerDataC3D = [-markerDataC3D(:, 1), markerDataC3D(:, 3), markerDataC3D(:, 2)]; %change reference frame
    
    markerDataC3D = markerDataC3D(1:min(length(markerDataC3D), length(markerDataCustomLog)), :);
    markerDataCustomLog = markerDataCustomLog(1:min(length(markerDataC3D), length(markerDataCustomLog)), :);
    
    %% pick out a good stretch of the first column of the custom log marker without nan's
    customLogDataToSyncWithTrackedC3D = markerDataCustomLog(:, 1);
    trackedDataToSyncWithCustomLog = markerDataC3D(:, 1);
    
    if (sum(isnan(customLogDataToSyncWithTrackedC3D))) % find a good stretch of data to use for sync check
        isGoodVector = ~isnan(customLogDataToSyncWithTrackedC3D);
        indecesOfBeginningsOfGoodStretches = find(diff(isGoodVector) > 0.5) + 1;
        indecesOfEndsOfGoodStretches = find(diff(isGoodVector) < -0.5);
        if (isempty(indecesOfEndsOfGoodStretches))
            indecesOfEndsOfGoodStretches = length(customLogDataToSyncWithTrackedC3D) - 1;
        end
        
        if (indecesOfEndsOfGoodStretches(1) < indecesOfBeginningsOfGoodStretches(1))
            indecesOfBeginningsOfGoodStretches = [1; indecesOfBeginningsOfGoodStretches];
        end
        
        if (indecesOfEndsOfGoodStretches(end) < indecesOfBeginningsOfGoodStretches(end))
            indecesOfEndsOfGoodStretches = ...
                [indecesOfEndsOfGoodStretches; length(customLogDataToSyncWithTrackedC3D)];
        end
        
        lengthsOfGoodStretches = indecesOfEndsOfGoodStretches - indecesOfBeginningsOfGoodStretches;
        [lengthOfLongestStretch, indexOfMaxStretch] = max(lengthsOfGoodStretches);
        indecesOfGoodData = (5:(lengthOfLongestStretch - 5)) + indecesOfBeginningsOfGoodStretches(indexOfMaxStretch); % take longest good stretch, minus 5 samples from each end
        if (sum(isnan(customLogDataToSyncWithTrackedC3D(indecesOfGoodData))))
            error('failed to find a period of custom log data without nans'); % something would be very wrong
        end
        
        customLogDataToSyncWithTrackedC3D = customLogDataToSyncWithTrackedC3D(indecesOfGoodData);
        trackedDataToSyncWithCustomLog = trackedDataToSyncWithCustomLog(indecesOfGoodData);
    end
    
    firstIndexToUseInSync = 1;
    indecesToUseForSync = firstIndexToUseInSync:length(trackedDataToSyncWithCustomLog);
    
    %% check sync/lag with xcorr
    C3DdataRange = max(max(markerDataC3D) - min(markerDataC3D));
    if (C3DdataRange < 0.1)
        if ignoreStillTrials
            if ~c3dUnitsAreInMeters
                disp(['Warning: this c3d (' c3dFilenameToLoad ') is either barely moving or in meters instead of mm and the c3dUnitsAreInMeters flag was not set. Processing is continuing due to ignoreStillTrials flag'])
            end
        else
            if ~c3dUnitsAreInMeters
                error(sprintf(['It looks like the marker data from the c3d is either barely moving, or it''s not in mm, it''s in meters. Check this... \n' ...
                    'If the c3d is actually in meters, rerun this code adding to the function call: savePhasespaceAndAnalogDataInC3D(..., ''c3dUnitsAreInMeters'', 1)\n' ...
                    'If the c3d is a still trial, rerun this code adding to the function call: savePhasespaceAndAnalogDataInC3D(..., ''ignoreStillTrials'', 1)']));
            else
                error(sprintf(['It looks like the marker date from the c3d is barely moving.'...
                    'If the c3d is a still trial, rerun this code adding to the function call: savePhasespaceAndAnalogDataInC3D(..., ''ignoreStillTrials'', 1)']));
            end
        end
    end
    
    % improper scaling will do weird things here.
    [vals, lags] = xcorr( ...
        (customLogDataToSyncWithTrackedC3D - mean(customLogDataToSyncWithTrackedC3D)), ...
        (trackedDataToSyncWithCustomLog(indecesToUseForSync) - mean(trackedDataToSyncWithCustomLog(indecesToUseForSync))), ...
        'none');
    
    [maxCorr, maxCorrInd] = max(vals);
    maxLag = lags(maxCorrInd);
    
    differencesBetweenC3DMarkersAndCustomLog = markerDataCustomLog - markerDataC3D;
    %   if (max(max(abs(differencesBetweenC3DMarkersAndCustomLog))) > 1e-2)
    %     fprintf('custom log and c3d log give different results,
    %     performing time lag rectification, calculated max lag is %g\n',
    %     maxLag); if (abs(maxLag) > 1000)
    %       error('the lag between the c3d marker data and the
    %       markerDataCustomLog is large, and the resulting curves don''t
    %       match well. Are you sure you are specifying c3d and customLogs
    %       that correspond to each other?');
    %     end
    %   end
    
    customLogIndecesToUse = 1:length(markerDataCustomLog);
    customLogIndecesToUse = customLogIndecesToUse - maxLag;
    
    customLogIndecesToUse(customLogIndecesToUse > length(markerDataCustomLog)) = [];
    customLogIndecesToUse(customLogIndecesToUse < 1) = [];
    
    %% optionally fill gaps in data
    if (shouldFillGapsInMarkerData)
        for markerNumber = 1:length(trackedC3dTrial.markerNames)
            markerName = trackedC3dTrial.markerNames{markerNumber};
            if (sum(sum(isnan(trackedC3dTrial.(markerName)))) > 0)
                fprintf('warning, marker %s isn''t fully tracked in the c3d %s\n', markerName, c3dFilenameToLoad);
                if (length(varargin) == 0)
                    [trackedData] = fillGapsInMarkerData(markerName, trackedC3dTrial.(markerName));
                else
                    [trackedData] = fillGapsInMarkerData(markerName, trackedC3dTrial.(markerName), varargin{:});
                end
                trackedC3dTrial.(markerName) = trackedData;
            end
        end
    end
    
    %% swap in the c3d data for the custom log data, optionally truncate
    badMarkerNumbers = [];
    for markerNumber = 1 : length(markerNames)
        %     markerNumber
        markerName = markerNames{markerNumber};
        if (strcmp(markerName, 'frameNumber'))
            continue;
        end
        
        if (c3dMarkerNamesHaveThreePaddedZeros)
            markerNameC3DCorrespondingToPhasespace= sprintf('M%03g', str2double(markerName(2:end)));
        else
            markerNameC3DCorrespondingToPhasespace = sprintf(markerName);
        end
        
        %     markerName
        
        try
            try
                markerDataC3D = trackedC3dTrial.(markerNameC3DCorrespondingToPhasespace);
            catch e
                markerDataC3D = trackedC3dTrial.(upper(markerNameC3DCorrespondingToPhasespace));
            end
            markerDataC3D = markerDataC3D(customLogIndecesToUse, :);
            
            markerDataC3D = [-markerDataC3D(:, 1), markerDataC3D(:, 3), markerDataC3D(:, 2)];
            
            if (isfield(customPhasespaceLog, markerName))
                markerDataCustomLog = customPhasespaceLog.(markerName);
                markerDataC3D = markerDataC3D(1:min(length(markerDataC3D), length(markerDataCustomLog)), :);
            end
            
            if (sum(sum(isnan(markerDataC3D))))
                disp([markerName ' has nans'])
                if (errorOnFindingNaNs)
                    error('shouldn''t replace marker data with nans!');
                end
            end
            
            % time correction
            markerLinSpace = 1:length(markerDataC3D);
            
            % never write a comment like the following. Where, for example,
            % might one find "his plots"? Does he intend to have these
            % available for the next few decades anytime anyone wants to
            % see them? And if someone comes up to him and asks him for
            % "his plots", would he know what to show them?
            actuallySyncedMarkerTimes = markerLinSpace * (1+(3.328e-5)) + (12.74/analogToMarkerRatio); % computed by Justin Sung, see his plots for justification
            actuallySyncedMarkerTimes = actuallySyncedMarkerTimes(actuallySyncedMarkerTimes < max(markerLinSpace));
            markerDataInterp = interp1(markerLinSpace, markerDataC3D, actuallySyncedMarkerTimes);
            
            if any(truncateFrames)>0
                if length(analogTruncate)>1
                    markerDataInterp = markerDataInterp(((analogTruncate(1)-1)/analogToMarkerRatio+1):(analogTruncate(2)/analogToMarkerRatio),:);
                else
                    markerDataInterp = markerDataInterp(((analogTruncate(1)-1)/analogToMarkerRatio+1):end,:);
                end
            end
            
            customPhasespaceLog.(markerName) = markerDataInterp;% * 1000.0;
            
        catch e
            %       e
            fprintf('had a problem trying to process marker: %s\n', markerName);
            badMarkerNumbers(length(badMarkerNumbers) + 1) = markerNumber;
        end
    end
    markerNames(badMarkerNumbers) = [];
    numberOfMarkers = length(markerNames);
end

% At this point the marker data may have a shorter length due to
% truncating datapoints that exceed original data from the sampling
% rate alignment adjustment. If so, we further cut both.

shortestLength = min([length(customPhasespaceLog.(markerNames{find(~strcmp(markerNames,'frameNumber'),1)})) * analogToMarkerRatio, ...
    length(analogChannelData{1})]);
shortestLength = shortestLength - mod(shortestLength, analogToMarkerRatio);

for i = 1:length(analogChannelData)
    analogChannelData{i} = ...
        analogChannelData{i}(1:shortestLength);
end

for markerNumber = 1:numberOfMarkers
    markerName = markerNames{markerNumber};
    if (strcmp(markerName, 'frameNumber'))
        continue;
    end
    markerData = customPhasespaceLog.(markerName);
    customPhasespaceLog.(markerName) = ...
        markerData(1 : shortestLength / analogToMarkerRatio, :);
end

if (diff([length(customPhasespaceLog.(markerName)) * analogToMarkerRatio, ...
        length(analogChannelData{1})]))
    error('trimmed analog data length doesn''t match trimmed marker data length again'); %shouldn't really happen
end

%% Load Perturbation File Data
if ~isempty(perturbationFile)
    %Handle File Types for perturbationFile
    
    if strcmp(perturbationFile(end-4:end),'.tdms')
        pertmatfilename = [perturbationFile(1:end-5) '.mat'];
        if exist(pertmatfilename,'file')
            %If perturbationFile has been processed before, look for .mat
            %file in same folder as the tdms file.  This is done because
            %reading .tdms files in matlab is very slow
            load(pertmatfilename,'P')
        else
            fprintf(1,'Reading .tdms perturber file\n');
            P = ReadPerturberData(perturbationFile);
            save(pertmatfilename,'P');
        end
        
        loopCount = P.LoopCount;
        validTimes=loopCount>0;
        t=loopCount(validTimes);
        
        Perturber1Force=P.Force1(validTimes);
        Perturber1Encoder=P.Encoder1(validTimes);
        Perturber1Command = P.Command1(validTimes);
        Perturber1SetPoint = P.SetPoint1(validTimes);
        
        Perturber2Force=P.Force2(validTimes);
        Perturber2Encoder=P.Encoder2(validTimes);
        Perturber2Command = P.Command2(validTimes);
        Perturber2SetPoint = P.SetPoint2(validTimes);
        
    elseif strcmp(perturbationFile(end-4:end),'.xlsx')
        fprintf(1,'Reading .xlsx perturber file \n');
        perturber1=xlsread(perturbationFile,'Perturber 1');
        perturber2=xlsread(perturbationFile,'Perturber 2');
        sysTime=xlsread(perturbationFile,'System');
        
        loopCount=sysTime(:,2);
        validTimes=loopCount>0;
        t=sysTime(validTimes,2);
        
        Perturber1Force=perturber1(validTimes,1);
        Perturber1Encoder=perturber1(validTimes,2);
        Perturber1Command = perturber1(validTimes,3);
        Perturber1SetPoint = perturber1(validTimes,4);
        
        Perturber2Force=perturber2(validTimes,1);
        Perturber2Encoder=perturber2(validTimes,2);
        Perturber2Command = perturber2(validTimes,3);
        Perturber2SetPoint = perturber2(validTimes,4);
    else
        warning(['Unknown file type for perturbation file: ' perturbationFile]);
        fprintf(1,'Not performing sync to perturbation file\n');
    end
    %% Sync Motion & Force Plate Data to Perturbation Data
    timecode = analogChannelData{14};
    
    ratioFor100 = 960/100;
    inds = find(~(timecode > 0.5));
    ind = find(diff(inds) >= 239, 1); % find a start point after an appropriate pause (low signal) in time code
    ia = inds(ind)+241;
    indicesForTimingSignalData = ia + round(ratioFor100 * ((30:6:96)-26)); % get 3 bytes of frame number in correct endian order
    frameNumberBinary = timecode(indicesForTimingSignalData);
    frameNumber = bin2dec(regexprep(num2str(frameNumberBinary > 0.5), ' ', ''));
    
    frameNumberBinaryNext = timecode(indicesForTimingSignalData + 960); % double check this
    frameNumberNext = bin2dec(regexprep(num2str(frameNumberBinaryNext > 0.5), ' ', ''));
    
    if frameNumberNext ~= frameNumber + 1
        error('Perturber sync problem')
    end
    
    % here frames is just the number of elapsed seconds on the perturber,
    % get index of equal number of whole seconds
    forceIndex = (inds(ind)+1):floor((length(timecode)-(inds(ind)+1))/960)*960+(inds(ind));
    secondsElapsedInForceFile = length(forceIndex)/960;
    secondsRemainingInPerturbFile = floor(length(t((frameNumber*100+1):end))/100);
    
    if secondsRemainingInPerturbFile >= secondsElapsedInForceFile
        % use force file end to truncate
        motionIndex  = ceil(forceIndex(1)/2):(ceil(forceIndex(end)/2)-1);
        perturbIndex = (frameNumber*100+1):((frameNumber+secondsElapsedInForceFile)*100);
    else
        forceIndex = (inds(ind)+1):((inds(ind)+960*secondsRemainingInPerturbFile));
        motionIndex  = ceil(forceIndex(1)/2):(ceil(forceIndex(end)/2)-1);
        perturbIndex = (frameNumber*100+1):((frameNumber+secondsRemainingInPerturbFile)*100);
    end
    
    oldNumOfChannels = length(analogChannelData);
    %Discard non-synced analog channel frames
    for i = 1:oldNumOfChannels
        analogChannelData{i} = analogChannelData{i}(forceIndex);
    end
    
    % Add new analog channels corresponding to force and encoder readings.

    analogChannelData{oldNumOfChannels+1} = Perturber1Force(perturbIndex);
    analogChannelData{oldNumOfChannels+2} = Perturber2Force(perturbIndex);
    analogChannelData{oldNumOfChannels+3} = Perturber1Encoder(perturbIndex);
    analogChannelData{oldNumOfChannels+4} = Perturber2Encoder(perturbIndex);
    analogChannelData{oldNumOfChannels+5} = Perturber1Command(perturbIndex);
    analogChannelData{oldNumOfChannels+6} = Perturber2Command(perturbIndex);
    analogChannelData{oldNumOfChannels+7} = Perturber1SetPoint(perturbIndex);
    analogChannelData{oldNumOfChannels+8} = Perturber2SetPoint(perturbIndex);
    
    %Add Names for these force channels
    analogChannelNames{oldNumOfChannels+1} = 'F1';
    analogChannelNames{oldNumOfChannels+2} = 'F2';
    analogChannelNames{oldNumOfChannels+3} = 'E1';
    analogChannelNames{oldNumOfChannels+4} = 'E2';
    analogChannelNames{oldNumOfChannels+5} = 'C1';
    analogChannelNames{oldNumOfChannels+6} = 'C2';
    analogChannelNames{oldNumOfChannels+7} = 'S1';
    analogChannelNames{oldNumOfChannels+8} = 'S2';
    
    %Discard non-synced motion data frames
    for i = 1:length(markerNames)
        customPhasespaceLog.(markerNames{i}) =  customPhasespaceLog.(markerNames{i})(motionIndex,:);
    end
end

%% Switch Force Plates
if SwitchForcePlates
   analogChannelData(1:12) = analogChannelData([7:12 1:6]); 
end
%% Filter Forces
if ForceFilter
    LeftChannels = [1 2 3 4 5 6];
    RightChannels = [7 8 9 10 11 12];
    Lgrf = horzcat(500*analogChannelData{LeftChannels(1)}', ...
        -500*analogChannelData{LeftChannels(2)}', ...
        1000*analogChannelData{LeftChannels(3)}');
    Rgrf = horzcat(500*analogChannelData{RightChannels(1)}', ...
        -500*analogChannelData{RightChannels(2)}', ...
        1000*analogChannelData{RightChannels(3)}');
    
    LMom = horzcat(820*analogChannelData{LeftChannels(4)}', ...
        260*analogChannelData{LeftChannels(5)}', ...
        400*analogChannelData{LeftChannels(6)}');
    
    
    RMom = horzcat(800*analogChannelData{RightChannels(4)}', ...
        250*analogChannelData{RightChannels(5)}', ...
        400*analogChannelData{RightChannels(6)}');
    
    [B,A] = butter(2,ForceFilterFreq/(960/2));
    Rgrf = filtfilt(B,A,Rgrf);
    Lgrf = filtfilt(B,A,Lgrf);
    RMom = filtfilt(B,A,RMom);
    LMom = filtfilt(B,A,LMom);
    
    analogChannelData{LeftChannels(1)} = Lgrf(:,1)/500;
    analogChannelData{LeftChannels(2)} = Lgrf(:,2)/-500;
    analogChannelData{LeftChannels(3)} = Lgrf(:,3)/1000;
    
    analogChannelData{RightChannels(1)} = Rgrf(:,1)/500;
    analogChannelData{RightChannels(2)} = Rgrf(:,2)/-500;
    analogChannelData{RightChannels(3)} = Rgrf(:,3)/1000;
    
    analogChannelData{LeftChannels(4)} = LMom(:,1)/820;
    analogChannelData{LeftChannels(5)} = LMom(:,2)/260;
    analogChannelData{LeftChannels(6)} = LMom(:,3)/400;
    
    analogChannelData{RightChannels(4)} = RMom(:,1)/800;
    analogChannelData{RightChannels(5)} = RMom(:,2)/250;
    analogChannelData{RightChannels(6)} = RMom(:,3)/400;
end

%% Correct Force Plate Drift
if ForcePlateCorrection
    
    if ~exist('Rgrf','var')
        LeftChannels = [1 2 3 4 5 6];
        RightChannels = [7 8 9 10 11 12];
        Lgrf = horzcat(500*analogChannelData{LeftChannels(1)}', ...
            500*analogChannelData{LeftChannels(2)}', ...
            1000*analogChannelData{LeftChannels(3)}');
        Rgrf = horzcat(500*analogChannelData{RightChannels(1)}', ...
            500*analogChannelData{RightChannels(2)}', ...
            1000*analogChannelData{RightChannels(3)}');
    end
   
    %Find Heel-Strikes & Toe-Offs
    if strcmp(activity,'Run')
        [HS,TO] = findEventIndicesForRunningOnTreadmill(Rgrf, Lgrf, 960);
        
    elseif strcmp(activity,'Walk')
        [RHS,LHS,RTO,LTO] = findEventIndicesForWalkingOnTreadmill(Rgrf, Lgrf, 960);
        HS = sort([RHS;LHS]); TO = sort([RTO;LTO]);
    end
   

    [Rgrf,Lgrf,POS,VEL,ACC,COMWR] = UndriftForcePlates(Rgrf,Lgrf,HS,TO,960,...
        'mass',mass,'PolyFitOrder',PolyFitOrder,'verbose',verbose);
    
    if 1
        % For each region demarcated by gait events, check whether signal
        % is close to zero, and if so, add it's signal to the other force
        % plate and then set it equal to zero.
        events = sort([HS;TO]);
        SmallTol = 20;
        
        for i = 1:length(events)-1
            dex = events(i):events(i+1);
            
            if mean(Rgrf(dex,3)) < SmallTol
                RightIsSmall = true;
            else
                RightIsSmall = false;
            end
            if mean(Lgrf(dex,3)) < SmallTol
                LeftIsSmall = true;
            else
                LeftIsSmall = false;
            end
            
            if RightIsSmall && LeftIsSmall
                Rgrf(dex,:) = 0;
                Lgrf(dex,:) = 0;
            elseif RightIsSmall && ~LeftIsSmall
                Lgrf(dex,:) = Lgrf(dex,:) + Rgrf(dex,:);
                Rgrf(dex,:) = 0;
            elseif ~RightIsSmall && LeftIsSmall
                Rgrf(dex,:) = Lgrf(dex,:) + Rgrf(dex,:);
                Lgrf(dex,:) = 0;
            end
        end
        
    end
    
    analogChannelData{LeftChannels(1)} = Lgrf(:,1)/500;
    analogChannelData{LeftChannels(2)} = Lgrf(:,2)/-500;
    analogChannelData{LeftChannels(3)} = Lgrf(:,3)/1000;
    
    analogChannelData{RightChannels(1)} = Rgrf(:,1)/500;
    analogChannelData{RightChannels(2)} = Rgrf(:,2)/-500;
    analogChannelData{RightChannels(3)} = Rgrf(:,3)/1000;
    
    %% Create new Analog Channels for Position, Velocity, and Acceleration
%     oldNumOfChannels = length(analogChannelData);
%     analogChannelData{oldNumOfChannels+1} = POS(:,1);
%     analogChannelData{oldNumOfChannels+2} = POS(:,2);
%     analogChannelData{oldNumOfChannels+3} = POS(:,3);
%     analogChannelData{oldNumOfChannels+4} = VEL(:,1);
%     analogChannelData{oldNumOfChannels+5} = VEL(:,2);
%     analogChannelData{oldNumOfChannels+6} = VEL(:,3);
%     analogChannelData{oldNumOfChannels+7} = ACC(:,1);
%     analogChannelData{oldNumOfChannels+8} = ACC(:,2);
%     analogChannelData{oldNumOfChannels+9} = ACC(:,3);
%     
%     %Add Names for these force channels
%     analogChannelNames{oldNumOfChannels+1} = 'POSX';
%     analogChannelNames{oldNumOfChannels+2} = 'POSY';
%     analogChannelNames{oldNumOfChannels+3} = 'POSZ';
%     analogChannelNames{oldNumOfChannels+4} = 'VELX';
%     analogChannelNames{oldNumOfChannels+5} = 'VELY';
%     analogChannelNames{oldNumOfChannels+6} = 'VELZ';
%     analogChannelNames{oldNumOfChannels+7} = 'ACCX';
%     analogChannelNames{oldNumOfChannels+8} = 'ACCY';
%     analogChannelNames{oldNumOfChannels+9} = 'ACCZ';
    
end

%% don't ever ever name variables like this, name them something sane.
%   vfrt     = video frame rate, must be > 0, default 120 nfram    = number
%   of frames, default 500 nmkr     = number of markers, not < 0, default
%   20 avr      = analog to video ratio, default 13 achn     = number of
%   analog channels, not < 0, default 28 ftype    = Intel(1), DEC(2),
%   SGI(3), default 1 dtype    = Integer(1), Floating Point(2), default 2
%   pscal    = point data scaling factor, must be > 0, default 0.1

c3dServer = actxserver('C3DServer.C3D');
name = outputFilename; %'output.c3d';
videoFrameRate = phasespaceSampleFrequency;
numberOfFrames = length(customPhasespaceLog.(markerNames{find(~strcmp(markerNames, 'frameNumber'),1)}));
analogToVideoRatio = analogToMarkerRatio;
numberOfAnalogChannels = length(analogChannelData);
filetype = 3; %1; %
dataType = 2; %
pointDataScaleFactor = 1.0;

createc3d(c3dServer, name, videoFrameRate, numberOfFrames, numberOfMarkers, ... - 1, ...
    analogToVideoRatio, numberOfAnalogChannels, filetype, dataType, pointDataScaleFactor);

%% set the length of the trial explicitly, so that we can have trials longer than 2^16 frames
% parameterGroup = 'TRIAL'; c3dServer.AddGroup (0, parameterGroup, '',
% '0');
%
% parameterName = 'ACTUAL_END_FIELD'; description = ''; lockParam = '0';
% dataType = int16(2); %2); %4; %-1; %-1 : Character, 1 : Byte, 2 :
% Integer, 4 : Floating point numberOfDataDimensions = int16(2); %int16(2);
% %2; %2; %needs to be 2 for a string, since first dimension holds the
% length of each string dataDimensions = int16([1 2]); %1; %{uint8(1)}; %[1
% 1]; %[1, numberOfMarkers]; % [lengthOfEachString, numberOfStrings] data =
% int16([mod(numberOfFrames, (2^16 - 1)), floor(numberOfFrames / (2^16 -
% 1))]);
%
% c3dServer.AddParameter(parameterName, description, ...
%   parameterGroup, lockParam, dataType, numberOfDataDimensions, ...
%   dataDimensions, data); %[]); %
%
% parameterName = 'ACTUAL_START_FIELD'; data = int16([1, 0]);
%
% c3dServer.AddParameter(parameterName, description, ...
%   parameterGroup, lockParam, dataType, numberOfDataDimensions, ...
%   dataDimensions, data); %[]); %
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
    
    %   frameNums = 1:(numberOfFrames); for frameNum = frameNums
    %     for subFrameNum = 1:analogToVideoRatio
    %       c3dServer.SetAnalogData(channelNumber, frameNum, subFrameNum,
    %       sin(frameNum/10));
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
    %   markerIndex markerName
    
    % this is how it used to be, but we probably instead want to use the
    % data from the c3d?
    %   markerData = customPhasespaceLog.(markerName);
    markerData = customPhasespaceLog.(markerName);
    
    if (sum(sum(isnan(markerData))))
        if (errorOnFindingNaNs)
            error('shoulnd''t output nan data from filled in c3d!');
        end
    end
    %   fprintf('correctly removed nans from marker %s\n', markerName);
    
    
    %   if (c3dUnitsAreInMeters)
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
            % we want to do a transformation: y->z -z->y x->x
            if (channelToWrite == 0)
                channelToWrite = 0;
            elseif (channelToWrite == 1)
                channelToWrite = 2;
            elseif (channelToWrite == 2)
                channelToWrite = 1;
                dataToWrite = -1.0 * dataToWrite;
            end
            
            %       if(fillMarkerGapsWithSplines)
            %         goodIndeces = ~isnan(dataToWrite); allIndeces =
            %         1:length(dataToWrite); dataToWrite =
            %         spline(allIndeces(goodIndeces),
            %         dataToWrite(goodIndeces), allIndeces);
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
        %       if (channel == 3) % residual channel... very important that
        %       this be non-negative
        %         c3dServer.SetPointData(markerNumber, channel, frameNum,
        %         0.0001);
        %       else % x,y,z channels
        %         c3dServer.SetPointData(markerNumber, channel, frameNum,
        %         sin(frameNum/10));
        %       end
        %     end
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


savec3d(c3dServer,name); %, 'output.c3d'); %
closec3d(c3dServer);

phasespaceAndTreadmillStructure.markers = customPhasespaceLog;

for i = 1:length(analogChannelNames)
    analogChannelName = analogChannelNames{i}; %(((end-5):end)); %
    phasespaceAndTreadmillStructure.analogChannels.(analogChannelName) = ...
        analogChannelData{i};
end

%%
a = loadPhasespaceRecord(outputFilename, 'forceReprocess', 1);

if ForcePlateCorrection
    NAnalog = numberOfFrames/phasespaceSampleFrequency*960;
    POS = POS(1:NAnalog,:);
    VEL = VEL(1:NAnalog,:);
    ACC = ACC(1:NAnalog,:);
    
    matfilename = createMatOutputFilename(outputFilename);
    savenames = {'POS' 'VEL' 'ACC' 'HS' 'TO'};
    save(matfilename,savenames{:},'-append');
end
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


