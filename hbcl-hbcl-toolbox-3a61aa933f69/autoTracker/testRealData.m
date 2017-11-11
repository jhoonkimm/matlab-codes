

%%
%
% test on a bit of real data
%

% clear all;
% close all;
format compact;

%%
mocapSkipFactor = 4;
mocapFrequency = 480;

componentsToPlot = 3; %1; %2; %

rightFootMarkerNames = {'M000', 'M001', 'M002', 'M003', 'M004'};
rightThighMarkerNames = {'M008', 'M009', 'M011', 'M029', 'M031'};

markerNames = rightThighMarkerNames; %rightFootMarkerNames; %

generateTemplateFromDataDuringTrial = 0; %1; %

directory = 'C:\\Users\\jrebula\\myProjects\\IMUGaitAnalysis\\rawData\\LateralWalkingStabilityMechanisms\\20120919-mikePilot\\';
filename = [directory 'Trial-02.c3d.mat'];
filledFilename = [directory 'Trial-02 filled.c3d.mat'];
standingFilename = [directory 'staticstanding-01 filled.c3d.mat'];
% filename = [directory 'Trial-02.c3d'];
% filledFilename = [directory 'Trial-02 filled.c3d'];
% standingFilename = [directory 'staticstanding-01 filled.c3d'];

templateFilename = [directory 'modelTemplate.mdh'];

%%


%%
loadPhasespace = @(varName, filename) evalin('base', [varName ' = load(''' filename '''); ' varName ' = ' varName '.phasespaceRecord;']);

loadPhasespace('rawData', filename);
loadPhasespace('filledData', filledFilename);
loadPhasespace('standingData', standingFilename);

indecesToLookAt = 360:mocapSkipFactor:1200; %1600; %5000; %360:400; %
% indecesToLookAt = 760:mocapSkipFactor:1200; %1600; %5000; %360:400; %
% startIndex = (760 + 20*4);
% endIndex = startIndex + (55*4);
% indecesToLookAt = startIndex : mocapSkipFactor : endIndex; %1600; %5000; %360:400; %

% indecesToLookAt = 1:4:(length(rawData.M000) / 10);


%%
% pick out some data to play with, arrange into cell arrays and such

standingMarkers = zeros(3, length(markerNames));
for markerNumber = 1:length(markerNames)
  markerName = markerNames{markerNumber};
  standingMarkers(:, markerNumber) = mean(standingData.(markerName))';
end
standingMarkers = standingMarkers - repmat(standingMarkers(:,1), [1, length(markerNames)]);

rawMarkers = cell(length(indecesToLookAt), 1);
filledMarkers = cell(length(indecesToLookAt), 1);
rawMarkerDataToPlot = zeros(3 * length(markerNames), length(indecesToLookAt));
filledMarkerDataToPlot = zeros(3 * length(markerNames), length(indecesToLookAt));
for i = 1:length(indecesToLookAt)
  index = indecesToLookAt(i);
  rawMarkers{i} = zeros(3, length(markerNames));
  for markerNumber = 1:length(markerNames)
    markerName = markerNames{markerNumber};
    rawMarkers{i}(:, markerNumber) = rawData.(markerName)(index, :)';
    filledMarkers{i}(:, markerNumber) = filledData.(markerName)(index, :)';
    
    rawMarkerDataToPlot(((markerNumber-1)*3 + 1) : (markerNumber*3), i) = rawData.(markerName)(index, :)';
    filledMarkerDataToPlot(((markerNumber-1)*3 + 1) : (markerNumber*3), i) = filledData.(markerName)(index, :)';
  end
end

rawMarkerDataToPlot = rawMarkerDataToPlot';
filledMarkerDataToPlot = filledMarkerDataToPlot';


cleanupPlotString = ['set(gca, ''XTick'', []); set(gca, ''YTick'', []); box off; ' ...
  'set(gca, ''color'', ''none''); set(gcf, ''color'', ''w''); set(gca, ''XTickLabel'', []); ' ...
  ...  'set(gca, ''xcolor'', get(gcf, ''color''));' ...
  'set(gca, ''xcolor'', get(gcf, ''color''));'
  ];
cleanupPlot = @() eval(cleanupPlotString);

signalsFig = figure;
for markerNumber = 1:length(markerNames)
  subplot(length(markerNames), 1, markerNumber );
  
  indecesToPlot = ((markerNumber-1)*3 + 1) : (markerNumber*3);
  indecesToPlot = indecesToPlot(componentsToPlot);
  
  plot(filledMarkerDataToPlot(:,indecesToPlot), 'k-', 'LineWidth', 4);
  hold on;
  plot(rawMarkerDataToPlot(:, indecesToPlot), 'r-', 'LineWidth', 4);
  cleanupPlot();
  
  ylabel(markerNames{markerNumber});
end

drawnow;

%%
%
%

r = RigidBody;
for i = 1:length(markerNames)
  r = r.addTrackerPositionInBodyFrame(standingMarkers(:,i) / 1000.0);
end

if (0)
  trajectoriesFig = figure;
  % subplot(2, 2, 1);
  for i = 1:1:length(indecesToLookAt)
    r.plotTheseWorldPoints(filledMarkers{i}, {'go'});
    r.plotTheseWorldPoints(rawMarkers{i}, {'ro'});
  end
  axis equal;
end

%%
position = rawMarkers{1}(:,1) / 1000.0;
orientation = [0 0 0]';
r.position = position;
r.orientation = orientation;
estimatedMarkerDataToPlot = rawMarkerDataToPlot * NaN;
% figure(trajectoriesFig);

printFrequency = 20;

estimatedPositions = zeros(3, length(indecesToLookAt));
estimatedOrientations = zeros(3, length(indecesToLookAt));

estimatedPositionsSlow = zeros(3, length(indecesToLookAt));
estimatedOrientationsSlow = zeros(3, length(indecesToLookAt));

observedTemplates = cell(1, length(indecesToLookAt));
templateFitErrors = zeros(1, length(indecesToLookAt));
tic
for i = 1 : length(indecesToLookAt)
  
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
  
  %   [positionSlow, orientationSlow, ...
  %     trackersFromTemplateSlow, templateFitError] = ...
  %     r.calculatePositionAndOrientationFromTrackersInWorld(rawMarkers{i} / 1000.0);
  
  [position, orientation, ...
    trackersFromTemplate, templateFitError] = ...
    r.calculatePositionAndOrientationFromTrackersInWorldLeastSquares(rawMarkers{i} / 1000.0);
  
  %   trackersFromTemplate = trackersFromTemplateSlow;
  
  %   max(max((trackersFromTemplateSlow - trackersFromTemplate) * 1000))
  
  r.position = position;
  r.orientation = orientation;
  
  %   position = position * 1000.0;
  trackersFromTemplate = trackersFromTemplate * 1000.0;
  
  templateFitErrors(i) = templateFitError;
  
  estimatedOrientations(:, i) = orientation;
  estimatedPositions(:, i) = position;
  
  estimatedOrientationsSlow(:, i) = orientationSlow;
  estimatedPositionsSlow(:, i) = positionSlow;
  
  if (mod(i, printFrequency) == 0)
    time = toc;
    fprintf('finished %d frames, average time of %g s/frame (current running time of %g s, %g%% done)\n', ...
      i, time / i, time, 100*(i / length(indecesToLookAt)));
  end
  
  %   for markerNumber = 1:length(markerNames)
  %     estimatedMarkerDataToPlot(i, ((markerNumber-1)*3 + 1) : (markerNumber*3)) = trackersFromTemplate(:, markerNumber);
  %   end
  
  calculatedTrackers = r.getTrackerPositionsInWorld() * 1000.0;
  
  for markerNumber = 1:length(markerNames)
    estimatedMarkerDataToPlot(i, ((markerNumber-1)*3 + 1) : (markerNumber*3)) = trackersFromTemplate(:, markerNumber);
    %     estimatedMarkerDataToPlot(i, ((markerNumber-1)*3 + 1) : (markerNumber*3)) = calculatedTrackers(:, markerNumber);
  end
  
  
  %   r.position = position;
  %   r.orientation = orientation;
end


%%
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



%% filter templates...

if (generateTemplateFromDataDuringTrial)
  figure()
  plot(templateFitErrors)
  
  filterCuttoffFrequency = 10; %1; % hz
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
  
  figure;
  plot(reshape(stackedTemplatesRaw(:, 1, :), [size(stackedTemplatesRaw, 1), size(stackedTemplatesRaw, 3)])', 'r')
  hold on;
  plot(reshape(stackedTemplatesFiltered(:, 1, :), [size(stackedTemplatesFiltered, 1), size(stackedTemplatesFiltered, 3)])', 'k')
end




%% try smoothing the individually calculated frames:

% r = RigidBody;
% for i = 1:length(markerNames)
%   r = r.addTrackerPositionInBodyFrame(standingMarkers(:,i) / 1000.0);
% end


filterCuttoffFrequency = 10; %1; % hz
Wn = filterCuttoffFrequency / ((mocapFrequency / mocapSkipFactor) / 2);
filterOrder = 2;
[bButterworthFilter, aButterworthFilter] = butter(filterOrder, Wn);

estimatedOrientationsFiltered = estimatedOrientations;
estimatedPositionsFiltered = estimatedPositions;

for k = 1:3
  estimatedOrientationsFiltered(k, :) = filtfilt(bButterworthFilter, aButterworthFilter, estimatedOrientations(k, :));
  estimatedPositionsFiltered(k, :) = filtfilt(bButterworthFilter, aButterworthFilter, estimatedPositions(k, :));
  %   estimatedOrientationsFiltered(k, :) = estimatedOrientations(k, :);
  %   estimatedPositionsFiltered(k, :) = estimatedPositions(k, :);
end


estimatedFilteredMarkerDataToPlot = estimatedMarkerDataToPlot*NaN;
for i = 1:length(indecesToLookAt)
  
  r.position = estimatedPositionsFiltered(:, i);
  r.orientation = estimatedOrientationsFiltered(:, i);
  
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
  
  trackersFromTemplate = r.getTrackerPositionsInWorld() * 1000.0;  
  
  for markerNumber = 1:length(markerNames)
    estimatedFilteredMarkerDataToPlot(i, ((markerNumber-1)*3 + 1) : (markerNumber*3)) = trackersFromTemplate(:, markerNumber);
  end
end


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
    
    %   for markerNumber = 1:length(markerNames)
    %     estimatedMarkerDataToPlot(i, ((markerNumber-1)*3 + 1) : (markerNumber*3)) = trackersFromTemplate(:, markerNumber);
    %   end
    %
    %   r.position = position;
    %   r.orientation = orientation;
    
    %   toc
  end
  
  
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


%%

hAll = findall(0);
for idx = 1 : length(hAll)
  try
    set(hAll(idx), 'fontsize', 12);
  catch
  end
end



