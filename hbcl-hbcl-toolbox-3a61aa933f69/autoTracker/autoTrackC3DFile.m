function [allInterpolatedMarkers, rigidBodyStates, rawData] = autoTrackC3DFile(filename, standingFilename, templateFilename, varargin)

% directory = 'C:\\Users\\jrebula\\myProjects\\IMUGaitAnalysis\\rawData\\LateralWalkingStabilityMechanisms\\20120919-mikePilot\\';
% filename = [directory 'Trial-02.c3d.mat'];
% filledFilename = [directory 'Trial-02 filled.c3d.mat'];
% standingFilename = [directory 'staticstanding-01 filled.c3d.mat'];

c3dUnitsAreInMeters = 1;
mocapSkipFactor = 1;
mocapFrequency = 480;
componentsToPlot = 3; %1; %2; %
shouldPlot = 0; %1; %
printFrequency = 50000;
markerRenaming = {};
templateFitErrorMaximum = 0.015; %0.025; %0.005; %0.01; %0.03;
maxPointIndexToLookAt = []; %1000
verbose = 0;
passUninterpolatedMarkers = true;
selectSegments = ''; %'RTA, RPV, RTH, LTH';...'RPV,RTH,LTH,RSK,LSK';

% if you don't have, or don't want to specify, a templateFilename, you can
% instead specify these two variables to define which markers are
% associated with which rigid bodies. For example, 
% segmentNames = {'rightFoot', 'leftFoot'}
% rigidBodySegments.rightFoot = {'M000', 'M001', 'M002', 'M004'};
% rigidBodySegments.leftFoot = {'M010', 'M101', 'M012', 'M041'};
segmentNames = [];
rigidBodySegments = [];

for i = 1 : 2 : length(varargin)
  option = varargin{i};
  val = varargin{i + 1};
  switch option
    case 'c3dUnitsAreInMeters'
      c3dUnitsAreInMeters = value;
    case 'mocapSkipFactor'
      mocapSkipFactor = val;
    case 'mocapFrequency'
      mocapFrequency = val;
    case 'componentsToPlot'
      componentsToPlot = val;
    case 'shouldPlot'
      shouldPlot = val;
    case 'printFrequency'
      printFrequency = val;
    case 'markerRenaming'
      markerRenaming = val;
    case 'templateFitErrorMaximum'
      templateFitErrorMaximum = val;
    case 'maxPointIndexToLookAt'
      maxPointIndexToLookAt = val;
    case 'rigidBodySegments'
      rigidBodySegments = val;
    case 'segmentNames'
      segmentNames = val;
    case 'verbose'
      verbose = val;
    case 'passUninterpolatedMarkers'
      passUninterpolatedMarkers = val;
    case 'selectSegments'
      selectSegments = val;      
  end
end

clear allInterpolatedMarkers;
overallTimer = tic;


%%
generateTemplateFromDataDuringTrial = 0; %1; %

forceReprocessC3Ds = 0;
rawData = loadPhasespaceRecord(filename, ...
  'c3dUnitsAreInMeters', c3dUnitsAreInMeters, 'forceReprocess', forceReprocessC3Ds);
% rawData;

[rawData] = renameMarkersInC3DData(rawData, markerRenaming);

standingData = loadPhasespaceRecord(standingFilename, ...
  'c3dUnitsAreInMeters', c3dUnitsAreInMeters, 'forceReprocess', forceReprocessC3Ds);

if passUninterpolatedMarkers
  allInterpolatedMarkers=rawData;
  allInterpolatedMarkers=rmfield(allInterpolatedMarkers,'markerNames');
end

%%
numberOfIndecesToLookAt = (size(rawData.(rawData.markerNames{1}), 1));
if (~isempty(maxPointIndexToLookAt))
  numberOfIndecesToLookAt = min(numberOfIndecesToLookAt, maxPointIndexToLookAt);
end
indecesToLookAt = 1:mocapSkipFactor:numberOfIndecesToLookAt;

if (isempty(rigidBodySegments) || isempty(segmentNames))
  [rigidBodySegments, segmentNames] = parseRigidBodiesFromVisual3DModelTemplateFile(templateFilename);
  
%   fprintf('autotrackC3dFile... hack!\n');
%   rigidBodySegments.RTA = rigidBodySegments.RTA(3:end);
  
  if ~isempty(selectSegments)
    segmentsToUse=textscan(selectSegments,'%s','Delimiter',',');
    bodySegmentNames=fieldnames(rigidBodySegments);
    keepNames=false(size(segmentNames));
    keepBodySegments=false(size(bodySegmentNames));
    for segInc = 1:length(segmentsToUse{1})
      matchSegment = segmentsToUse{1}{segInc};
      keepNames=keepNames|strcmpi(segmentNames,matchSegment);
      keepBodySegments=keepBodySegments|strcmpi(bodySegmentNames,matchSegment); 
    end
    segmentNames=segmentNames(keepNames);
    rigidBodySegments=rmfield(rigidBodySegments,bodySegmentNames(~keepBodySegments));
  end
end

% fprintf('autotrackC3dFile: hack!!!!\n');
% segmentNames = ['RTA', segmentNames];
% segmentNames(end) = [];

cleanupPlotString = ['set(gca, ''XTick'', []); set(gca, ''YTick'', []); box off; ' ...
  'set(gca, ''color'', ''none''); set(gcf, ''color'', ''w''); set(gca, ''XTickLabel'', []); ' ...
  ...  'set(gca, ''xcolor'', get(gcf, ''color''));' ...
  'set(gca, ''xcolor'', get(gcf, ''color''));'
  ];
cleanupPlot = @() eval(cleanupPlotString);

for rigidBodyNumber = 1:length(segmentNames)
  segmentTimer = tic;
  timer = tic;
  markerNames = rigidBodySegments.(segmentNames{rigidBodyNumber});
  fprintf('processing rigid body %s\n', segmentNames{rigidBodyNumber});
  fprintf('markers on this body:\n');
  markerNames
  
  % average out standing marker locations for template body
  standingMarkers = zeros(3, length(markerNames));
  badMarkers = {};
  badMarkerIndeces = [];
  
  allMarkersKnownOriginal = fieldnames(standingData);
  allMarkersKnown = upper(allMarkersKnownOriginal);
  
  for markerNumber = 1:length(markerNames) % markers not part of a rigid body are thrown away
    markerName = markerNames{markerNumber};
    %     if (isfield(standingData, markerName))
    if (ismember(upper(markerName), allMarkersKnown))
      memberIndex = find(strcmp(upper(markerName), allMarkersKnown));
      markerName = allMarkersKnownOriginal{memberIndex};
      standingMarkers(:, markerNumber) = nanmean(standingData.(markerName))';
    else
      badMarkers{length(badMarkers) + 1} = markerName;
      badMarkerIndeces(length(badMarkerIndeces) + 1) = markerNumber;
    end
  end
  standingMarkers(:, badMarkerIndeces) = [];
  markerNames(badMarkerIndeces) = [];
  standingMarkers = standingMarkers - repmat(standingMarkers(:,1), [1, length(markerNames)]);
  
  %   fprintf('standing markers: \n');
  %   markerNames
  
  % build up block of trial marker data, filling in NaNs for nonexistent
  % markers
  
  allMarkersKnownOriginal = fieldnames(rawData);
  allMarkersKnown = upper(allMarkersKnownOriginal);
  
  rawMarkers = cell(length(indecesToLookAt), 1);
  allMarkersInBlock = zeros(length(indecesToLookAt), length(markerNames) * 3);
  for markerNumber = 1:length(markerNames)
    markerName = markerNames{markerNumber};
    %     if (isfield(rawData, markerNames{markerNumber}))
    if (ismember(upper(markerName), allMarkersKnown))
      memberIndex = find(strcmp(upper(markerName), allMarkersKnown));
      markerName = allMarkersKnownOriginal{memberIndex};
      
      allMarkersInBlock(:, ((markerNumber - 1) * 3) + (1:3)) = ...
        rawData.(markerName)(1:length(indecesToLookAt), :);
    else
      allMarkersInBlock(:, ((markerNumber - 1) * 3) + (1:3)) = ...
        zeros(length(indecesToLookAt), 3) * NaN;
    end
  end
  
  for i = 1:length(indecesToLookAt)
    index = indecesToLookAt(i);
    rawMarkers{i} = reshape(allMarkersInBlock(index, :), [3, length(markerNames)]);
  end
  
  rawMarkerDataToPlot = allMarkersInBlock; %rawMarkerDataToPlot';
  
  if (shouldPlot)
    signalsFig = figure;
    for markerNumber = 1:length(markerNames)
      subplot(length(markerNames), 1, markerNumber );
      
      indecesToPlot = ((markerNumber-1)*3 + 1) : (markerNumber*3);
      indecesToPlot = indecesToPlot(componentsToPlot);
      
      %     plot(filledMarkerDataToPlot(:,indecesToPlot), 'k-', 'LineWidth', 4);
      %     hold on;
      plot(rawMarkerDataToPlot(:, indecesToPlot), 'r-', 'LineWidth', 4);
      hold on;
      cleanupPlot();
      ylabel(markerNames{markerNumber});
    end
    drawnow;
  end
  
  r = RigidBody;
  for i = 1:length(markerNames)
    r = r.addTrackerPositionInBodyFrame(standingMarkers(:,i) / 1000.0);
  end
  
  %     trajectoriesFig = figure;
  position = rawMarkers{1}(:,1) / 1000.0;
  orientation = [0 0 0]';
  r.position = position;
  r.orientation = orientation;
  estimatedMarkerDataToPlot = rawMarkerDataToPlot * NaN;
  % figure(trajectoriesFig);
  
  estimatedPositions = zeros(3, length(indecesToLookAt));
  estimatedOrientations = zeros(3, length(indecesToLookAt));
%   estimatedSubsetOfTemplateMarkersForTransform = zeros(length(markerNames), length(indecesToLookAt));
  
  observedTemplates = cell(1, length(indecesToLookAt));
  templateFitErrors = zeros(1, length(indecesToLookAt));
  tic
  for i = 1:length(indecesToLookAt)
    if (generateTemplateFromDataDuringTrial)
      if (sum(sum(isnan(rawMarkers{i}))) == 0)
        %     if (sum(sum(isnan(rawMarkers{i})) < 0.5) >= 3)
        r = RigidBody;
        r.position = position; % / 1000.0;
        r.orientation = orientation;
        for j = 1:length(markerNames)
          r = r.addTrackerPositionInWorldFrame(rawMarkers{i}(:,j) / 1000.0);
        end
        observedTemplates{i} = r.trackersInBodyFrame;
      end
    else
      %       observedTemplates{i} = r.trackersInBodyFrame * NaN;
    end
    
    %       [position, orientation, ...
    %         trackersFromTemplate, templateFitError] = ...
    %         r.calculatePositionAndOrientationFromTrackersInWorld(rawMarkers{i} / 1000.0);
    %       try
    
    
    %%
    
%     allData = [];
%     for tickNumber = 1 : length(rawMarkers)
%     allData = [allData; rawMarkers{tickNumber}(1, :)];
%     end
    
    
    %%

    if (1 && sum(~isnan(rawMarkers{i}(1,:))) >= 3)
      [position, orientation, ...
        trackersFromTemplate, templateFitError] = ...
        r.calculatePositionAndOrientationFromTrackersInWorldLeastSquares(rawMarkers{i} / 1000.0, ...
        'shouldPlotFit', 0);
      
      if (templateFitError > templateFitErrorMaximum)
        fprintf('large error, knockout fitting, error = %g\n', templateFitError);
        goodMarkers = find(~isnan(rawMarkers{i}(1,:)));
        
        numberOfRawMarkers = length(goodMarkers); %size(rawMarkers{i}, 2);
        residualErrors = zeros(numberOfRawMarkers, 1);
        for markerToKnockout = 1:numberOfRawMarkers
          rawMarkersToUse = rawMarkers{i};
          rawMarkersToUse(:, goodMarkers(markerToKnockout)) = NaN;
          residualErrors(markerToKnockout) = NaN;
          if (sum(~isnan(rawMarkersToUse(1,:))) >= 3)
            [position, orientation, ...
              trackersFromTemplate, templateFitError] = ...
              r.calculatePositionAndOrientationFromTrackersInWorldLeastSquares(rawMarkersToUse / 1000.0);
            residualErrors(markerToKnockout) = templateFitError;
          end
        end
        
        [bestError, bestKnockoutMarker] = min(residualErrors);
        
        if ((bestError > templateFitErrorMaximum) || isnan(bestError))
          position = ones(3, 1) * NaN;
          orientation = ones(3, 1) * NaN;
          trackersFromTemplate = rawMarkers{i} * NaN;
          templateFitError = NaN;
          if (verbose > 1)
            fprintf('insufficient markers for fit, NaNing out tick %g!\n', i);
          end
        else
          rawMarkersToUse = rawMarkers{i};
          rawMarkersToUse(:, goodMarkers(bestKnockoutMarker)) = NaN;
          [position, orientation, ...
            trackersFromTemplate, templateFitError] = ...
            r.calculatePositionAndOrientationFromTrackersInWorldLeastSquares(rawMarkersToUse / 1000.0);
          if (verbose > 1)
            fprintf('found good knockout set, %g!\n', templateFitError);
          end
        end
      end
      
    else
      position = ones(3, 1) * NaN;
      orientation = ones(3, 1) * NaN;
      trackersFromTemplate = rawMarkers{i} * NaN;
      templateFitError = NaN;
      
      %         position = ones(3, 1);
      %         orientation = ones(3, 1);
      %         trackersFromTemplate = rawMarkers{i};
      %         templateFitError = 1;
    end
    %       catch e
    %         e
    %       end
    %   position = position * 1000.0;
    trackersFromTemplate = trackersFromTemplate * 1000.0;
    
    
    if (templateFitError > templateFitErrorMaximum)
      fprintf('large template fit error! This is a programming error if we get to this point!\n');
    end
    
    templateFitErrors(i) = templateFitError;
    
    estimatedOrientations(:, i) = orientation;
    estimatedPositions(:, i) = position;
    %     estimatedSubsetOfTemplateMarkersForTransform(:, i) = subsetOfTemplateMarkersForTransform;
    
    if (mod(i, printFrequency) == 0)
      time = toc(timer);
      fprintf('fitting rigid body poses: finished %d frames, average time of %g s/frame (current running time of %g s, %2.0f%% done)\n', ...
        i, time / i, time, 100*(i / length(indecesToLookAt)));
    end
    
    if (shouldPlot)
      for markerNumber = 1:length(markerNames)
        estimatedMarkerDataToPlot(i, ((markerNumber-1)*3 + 1) : (markerNumber*3)) = trackersFromTemplate(:, markerNumber);
      end
    end
    
    r.position = position;
    r.orientation = orientation;
  end
  
  % here we should really look at the fit errors... if some are really bad, then we should not trust it. We could 
  % rerun the fit procedure minus some markers whom we distrust, or we
  % could just throw out the bad fits and spline through them...
%   figure();
%   plot(templateFitErrors, 'k')
%   hold on;
%   plot(find(isnan(templateFitErrors)), 0.05 * ones(size(find(isnan(templateFitErrors)))), 'r*')
%   legend('errors', 'couldn''t find acceptable fit')
%   title(['templateFitErrors ' segmentNames{rigidBodyNumber}]);
    
  %
  allIndecesToInterp = 1:size(estimatedPositions, 2);
  nanIndecesToInterp = find(isnan(estimatedPositions(1,:)));
  nonNanIndecesToInterp = allIndecesToInterp;
  nonNanIndecesToInterp(nanIndecesToInterp) = [];
  
  try
    interpolatedEstimatedPostions = ...
      interp1(nonNanIndecesToInterp, estimatedPositions(:, nonNanIndecesToInterp)', allIndecesToInterp)';
  catch e
    e
  end
  interpolatedEstimatedOrientations = ...
    interp1(nonNanIndecesToInterp, estimatedOrientations(:, nonNanIndecesToInterp)', allIndecesToInterp)';
  
  estimatedPositions = interpolatedEstimatedPostions;
  estimatedOrientations = interpolatedEstimatedOrientations;
  
  %   if (sum(sum(isnan(estimatedPositions))))
  %     error('found some nan positions after interpolation!');
  %   end
  
  %
  if (shouldPlot)
    figure(signalsFig);
    for markerNumber = 1:length(markerNames)
      subplot(length(markerNames), 1, markerNumber);
      indecesToPlot = ((markerNumber-1)*3 + 1) : (markerNumber*3);
      indecesToPlot = indecesToPlot(componentsToPlot);
      plot(estimatedMarkerDataToPlot(:, indecesToPlot), 'g-', 'LineWidth', 2);
      hold on;
      cleanupPlot();
      ylabel(markerNames{markerNumber});
    end
    set(gca, 'xcolor', 'k');
    xlabel('time');
  end
  
  
  %% filter templates...
  if (generateTemplateFromDataDuringTrial)
    
    filterCuttoffFrequency = 230; %10; %1; % hz
    Wn = filterCuttoffFrequency / ((mocapFrequency / mocapSkipFactor) / 2);
    filterOrder = 2;
    [bButterworthFilter, aButterworthFilter] = butter(filterOrder, Wn);
    
    stackedTemplatesRaw = zeros( ...
      size(observedTemplates{1}, 1), ...
      size(observedTemplates{1}, 2), ...
      length(observedTemplates));
    stackedTemplatesFiltered = stackedTemplatesRaw;
    
    for i = 1:size(observedTemplates{1}, 1)
      for j = 1:size(observedTemplates{1}, 2)
        for k = 1:length(observedTemplates)
          stackedTemplatesRaw(i, j, k) = observedTemplates{k}(i, j);
        end
        stackedTemplatesFiltered(i,j,:) = filtfilt(bButterworthFilter, aButterworthFilter, stackedTemplatesRaw(i,j,:));
      end
    end
    
    if (shouldPlot)
      figure;
      plot(reshape(stackedTemplatesRaw(:, 1, :), [size(stackedTemplatesRaw, 1), size(stackedTemplatesRaw, 3)])', 'r')
      hold on;
      plot(reshape(stackedTemplatesFiltered(:, 1, :), [size(stackedTemplatesFiltered, 1), size(stackedTemplatesFiltered, 3)])', 'k')
    end
  end
  
  
  %% try smoothing the individually calculated frames:
  % r = RigidBody;
  % for i = 1:length(markerNames)
  %   r = r.addTrackerPositionInBodyFrame(standingMarkers(:,i) / 1000.0);
  % end
  
  filterCuttoffFrequency = 20; %230; %1; % hz
  Wn = filterCuttoffFrequency / ((mocapFrequency / mocapSkipFactor) / 2);
  filterOrder = 2;
  [bButterworthFilter, aButterworthFilter] = butter(filterOrder, Wn);
  
  estimatedOrientationsFiltered = estimatedOrientations;
  estimatedPositionsFiltered = estimatedPositions;
  
  for k = 1:3
    goodIndeces = ~isnan(estimatedOrientations(k, :));
    estimatedOrientationsFiltered(k, goodIndeces) = filtfilt(bButterworthFilter, aButterworthFilter, estimatedOrientations(k, goodIndeces));
    estimatedPositionsFiltered(k, goodIndeces) = filtfilt(bButterworthFilter, aButterworthFilter, estimatedPositions(k, goodIndeces));
    %       estimatedOrientationsFiltered(k, :) = estimatedOrientations(k, :);
    %       estimatedPositionsFiltered(k, :) = estimatedPositions(k, :);
  end
  
  for markerNumber = 1:length(markerNames)
    markerName = markerNames{markerNumber};
    allInterpolatedMarkers.(markerName) = zeros(length(indecesToLookAt), 3);
  end
  
  estimatedFilteredMarkerDataToPlot = estimatedMarkerDataToPlot * NaN;
  rigidBodyStates.(segmentNames{rigidBodyNumber}) = cell(length(indecesToLookAt), 1);
  tic
  for i = 1:length(indecesToLookAt)
    
    if (mod(i, printFrequency) == 0)
      time = toc;
      fprintf('generating marker positions from filtered rigid body poses: finished %d frames, average time of %g s/frame (current running time of %g s, %2.0f%% done)\n', ...
        i, time / i, time, 100*(i / length(indecesToLookAt)));
    end
    
    r.position = estimatedPositionsFiltered(:, i);
    r.orientation = estimatedOrientationsFiltered(:, i);
    
%     rigidBodyStates.(segmentNames{rigidBodyNumber}){i} = r;
    rigidBodyStates.(segmentNames{rigidBodyNumber}){i}.position = r.position;
    rigidBodyStates.(segmentNames{rigidBodyNumber}){i}.rotationMatrix = rodrigues(r.orientation);
    
    %   if (generateTemplateFromDataDuringTrial)
    %     if (sum(sum(isnan(rawMarkers{i}))) == 0)
    %       r = RigidBody;
    %       %       r.position = position / 1000.0;
    %       %       r.orientation = orientation;
    %       r.position = estimatedPositionsFiltered(:, i);
    %       r.orientation = estimatedOrientationsFiltered(:, i);
    %       for j = 1:length(markerNames)
    %         r = r.addTrackerPositionInWorldFrame(rawMarkers{i}(:,j) / 1000.0);
    %       end
    %     end
    %   end
    
    if (generateTemplateFromDataDuringTrial)
      r = RigidBody;
      r.position = estimatedPositionsFiltered(:, i);
      r.orientation = estimatedOrientationsFiltered(:, i);
      
      for j = 1:length(markerNames)
        %     r = r.addTrackerPositionInBodyFrame(observedTemplates{i}(:,j));
        r = r.addTrackerPositionInBodyFrame(stackedTemplatesFiltered(:,j,i));
      end
    end
    
    interpolatedTrackersFromTemplate = r.getTrackerPositionsInWorld() * 1000.0;
    
    for markerNumber = 1:length(markerNames)
      markerName = markerNames{markerNumber};
      %         allInterpolatedMarkers.(markerName)(:, i) = interpolatedTrackersFromTemplate(:, markerNumber);
      
%       fprintf('autotrackc3dFile: setting interpolated marker: %s\n', markerName);
      
      allInterpolatedMarkers.(markerName)(i, :) = interpolatedTrackersFromTemplate(:, markerNumber)';
      if (shouldPlot)
        estimatedFilteredMarkerDataToPlot(i, ((markerNumber-1)*3 + 1) : (markerNumber*3)) = interpolatedTrackersFromTemplate(:, markerNumber);
      end
    end
  end
  
  if (shouldPlot)
    figure(signalsFig);
    for markerNumber = 1:length(markerNames)
      subplot(length(markerNames), 1, markerNumber);
      indecesToPlot = ((markerNumber-1)*3 + 1) : (markerNumber*3);
      indecesToPlot = indecesToPlot(componentsToPlot);
      plot(estimatedFilteredMarkerDataToPlot(:, indecesToPlot), 'b-', 'LineWidth', 2);
      hold on;
      cleanupPlot();
      ylabel(markerNames{markerNumber});
    end
    set(gca, 'xcolor', 'k');
    xlabel('time');
    
    figure()
    plot(estimatedPositions');
    hold on;
    plot(estimatedPositionsFiltered', '--');
  end
  
  
  %% now try the smoothing over multiple gaps
  smoothOverGapsWithOptimization = 0;
  if (smoothOverGapsWithOptimization)
    
    smootherMovingWindowSize = 55;
    numberOfTicksToAddToSmootherPerIteration = 1; %4; %1;
    
    wholeStateEstimate = zeros(6, length(indecesToLookAt));
    
    currentMovingWindowStateEstimate = ...
      [estimatedPositions(:, 1:smootherMovingWindowSize); ...
      estimatedOrientations(:, 1:smootherMovingWindowSize)];
    
    estimatedInterpolatedMarkerDataToPlot = rawMarkerDataToPlot * NaN;
    
    printFrequency = 2;
    tic
    for i = (smootherMovingWindowSize + numberOfTicksToAddToSmootherPerIteration): ...
        numberOfTicksToAddToSmootherPerIteration : ...
        length(indecesToLookAt)
      
      indecesIntoAllData = ((i - smootherMovingWindowSize)+1) : i;
      
      if (generateTemplateFromDataDuringTrial)
        if (sum(sum(isnan(rawMarkers{indecesIntoAllData(1)}))) == 0)
          r = RigidBody;
          r.position = wholeStateEstimate(1:3, indecesIntoAllData(1)) / 1000.0;
          r.orientation = wholeStateEstimate(4:6, indecesIntoAllData(1));
          for j = 1:length(markerNames)
            r = r.addTrackerPositionInWorldFrame(rawMarkers{indecesIntoAllData(1)}(:,j) / 1000.0);
          end
        end
      end
      
      thisWindowMarkerTrajectories = cell(smootherMovingWindowSize, 1);
      for j = 1:length(thisWindowMarkerTrajectories)
        thisWindowMarkerTrajectories{j} = rawMarkers{indecesIntoAllData(j)} / 1000.0;
      end
      
      [interpolatedPositions, interpolatedOrientations, interpolatedMarkersEstimatedfromTemplate] = ...
        r.interpolateRigidBodyFromMarkers(thisWindowMarkerTrajectories);
      interpolatedPositions = interpolatedPositions * 1000;
      %   interpolatedMarkersEstimatedfromTemplate = interpolatedMarkersEstimatedfromTemplate * 1000.0;
      
      wholeStateEstimate(:, indecesIntoAllData) = [interpolatedPositions; interpolatedOrientations];
      
      for indexIntoAllDataNum = 1:length(indecesIntoAllData)
        indexIntoAllData = indecesIntoAllData(indexIntoAllDataNum);
        for markerNumber = 1:length(markerNames)
          estimatedInterpolatedMarkerDataToPlot(indexIntoAllData, ((markerNumber-1)*3 + 1) : (markerNumber*3)) = ...
            interpolatedMarkersEstimatedfromTemplate{indexIntoAllDataNum}(:, markerNumber);
        end
      end
      
      if (mod(i, printFrequency) == 0)
        time = toc;
        fprintf('finished %d frames, average time of %g s/frame (current running time of %g s, %g%% done)\n', ...
          i, time / i, time, 100*(i / length(indecesToLookAt)));
      end
    end
    
    if (shouldPlot)
      figure(signalsFig);
      for markerNumber = 1:length(markerNames)
        subplot(length(markerNames), 1, markerNumber);
        indecesToPlot = ((markerNumber-1)*3 + 1) : (markerNumber*3);
        indecesToPlot = indecesToPlot(componentsToPlot);
        plot(estimatedInterpolatedMarkerDataToPlot(:, indecesToPlot) * 1000, 'm-', 'LineWidth', 2);
        %   plot(estimatedInterpolatedMarkerDataToPlot(:, indecesToPlot), 'm-', 'LineWidth', 2);
        hold on;
        cleanupPlot();
        ylabel(markerNames{markerNumber});
      end
      set(gca, 'xcolor', 'k');
      xlabel('time');
    end
  end
  
  %%
  if (shouldPlot)
    hAll = findall(0);
    for idx = 1 : length(hAll)
      try
        set(hAll(idx), 'fontsize', 12);
      catch
      end
    end
  end
    fprintf('finished processing rigid body %s in %gs\n', segmentNames{rigidBodyNumber}, toc(segmentTimer));
end

totalTime = toc(overallTimer);
fprintf('finished processing autotracking in %gs, file: %s\n', totalTime, filename);



