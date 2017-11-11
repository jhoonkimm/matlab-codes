
close all;
% clear all;
clc;

%%
demoMode = 1; %0; %
forceReprocessEverything = 1; %0; %

%%
dataDir = 'C:\Users\jrebula\myProjects\IMUGaitAnalysis\rawData\overgroundVariabilityExperiment\20110512-PS and Cart Data_larkinsh\calibration\';

alignmentC3D = [dataDir '20110512_calibration-01.c3d'];
alignmentRST = [dataDir 'cartData_20110512_calibration1.rst'];

forceReprocessPhasespace = 0 || forceReprocessEverything;
alignmentPhasespaceRecord = loadPhasespaceRecord(alignmentC3D, 'verbose', 0, 'forceReprocess', forceReprocessPhasespace);
visualizePhasespaceRecord(alignmentPhasespaceRecord, 'plotInterval', 5, 'recordName', 'alignment trial');

alignmentMobileCartRecord = loadMobileCartRecord(alignmentRST, 'verbose', 0, 'forceReprocess', forceReprocessEverything);
[alignmentMobileCartMotion] = integrateMobileCartMotion(alignmentMobileCartRecord);

phasespaceCalibrationIndeces = []; %1:length(alignmentPhasespaceRecord.(alignmentPhasespaceRecord.markerNames{1}))/2;
[cartMocapAlignment] = getCartPhasespaceFrameAlignment(alignmentC3D, alignmentRST, ...
  'verbose', 20, 'forceReprocess', forceReprocessEverything, 'phasespaceCalibrationIndeces', phasespaceCalibrationIndeces);

[transformedAlignmentPhasespaceRecord] = transformPhasespaceRecordIntoWorldWithMobileCartAlignment(...
  alignmentPhasespaceRecord, alignmentMobileCartRecord, cartMocapAlignment);

figure;
visualizePhasespaceRecord(transformedAlignmentPhasespaceRecord, 'recordName', 'long walk trial', 'plotDropoutRateHistogram', 0);
visualizeMobileCartMotion(alignmentMobileCartMotion, 'recordName', 'alignment trial')
axis equal;
figure;
% makeVideoOfPhasespaceRecord('test2.avi', transformedAlignmentPhasespaceRecord, ...
%   'mobileCartMotion', alignmentMobileCartMotion);

% return;



%%
% dataDir = 'C:\Users\jrebula\myProjects\IMUGaitAnalysis\rawData\BLAST\20110819-MRL-I-S01\';
dataDir = 'C:\Users\jrebula\myProjects\IMUGaitAnalysis\rawData\BLAST\20110819-MRL-I-S01\';


% c3d files are from the phasespace mocap system
% rst files are from the mobile cart (encoders and integrated gyro)

alignmentC3D = [dataDir 'Cart-Alignment-02.c3d'];
alignmentRST = [dataDir 'Blast-cartAlignment1.rst'];

longWalkC3D = [dataDir 'longWalk\blastLongWalk-01.c3d'];
longWalkRST = [dataDir 'longWalk\blastLongWalk1.rst'];

% some things we want to be able to do:

% load in a c3d by itself
% look at the mocap record
% load in a cart file
% process cart record into cart motion
% look at the cart motion

% use a mocap record and cart record to calculate an alignment between the
% phasespace and cart frames

% transform mocap and cart record to world frame using an alignment
% look at mocap markers in world
% look at cart and mocap markers in world
% make video of mocap markers and cart in world, stationary viewpoint
% video, track target, dolly camera


%% The interface with the toolbox:


%% Phasespace record loading and visualization
% load in a c3d by itself
forceReprocessPhasespace = 0 || forceReprocessEverything;
alignmentPhasespaceRecord = loadPhasespaceRecord(alignmentC3D, 'verbose', 0, 'forceReprocess', forceReprocessPhasespace);
longWalkPhasespaceRecord = loadPhasespaceRecord(longWalkC3D, 'forceReprocess', forceReprocessPhasespace);

% look at the mocap record
if (demoMode)
  visualizePhasespaceRecord(alignmentPhasespaceRecord, 'plotInterval', 5, 'recordName', 'alignment trial');
  figure;
  visualizePhasespaceRecord(alignmentPhasespaceRecord, 'plotInterval', 5, 'recordName', 'long walk trial', ...
    'plotDropoutRateHistogram', 0);
  figure;
  numCols = 8;
  visualizePhasespaceRecord(alignmentPhasespaceRecord, 'plotInterval', 5, 'recordName', 'alignment trial', ...
    'axisHandlesToPlotOn', [subplot(2, numCols, [1:(numCols-1)]), subplot(2, numCols, numCols)]);
  visualizePhasespaceRecord(longWalkPhasespaceRecord, 'plotInterval', 50, 'recordName', 'long walk trial', ...
    'axisHandlesToPlotOn', [subplot(2, numCols, [(numCols+1):(2*numCols-1)]), subplot(2, numCols, 2*numCols)]);
end

%% Mobile cart record loading and visualization
% load in a cart file
alignmentMobileCartRecord = loadMobileCartRecord(alignmentRST, 'verbose', 0, 'forceReprocess', forceReprocessEverything);
longWalkMobileCartRecord = loadMobileCartRecord(longWalkRST, 'forceReprocess', forceReprocessEverything);

% process cart record into cart motion
[alignmentMobileCartMotion] = integrateMobileCartMotion(alignmentMobileCartRecord);
[longWalkMobileCartMotion] = integrateMobileCartMotion(longWalkMobileCartRecord);

% look at the cart motion
if (demoMode)
  figure;
  subplot(2, 1, 1);
  visualizeMobileCartMotion(alignmentMobileCartMotion, 'recordName', 'alignment trial')
  subplot(2, 1, 2);
  visualizeMobileCartMotion(longWalkMobileCartMotion, 'recordName', 'long walk trial', ...
    'orientationPlotInterval', 100)
end

%% Cart frame calibration
% use a mocap record and cart record to calculate an alignment between the
% phasespace and cart frames
% this calculates the alignment from loaded data:
% [cartMocapAlignment] = getCartPhasespaceFrameAlignment(alignmentPhasespaceRecord, alignmentMobileCartMotion);
% this calculates the alignment from filenames, and saves the results to a
% mat:
phasespaceCalibrationIndeces = []; %1:length(alignmentPhasespaceRecord.(alignmentPhasespaceRecord.markerNames{1}))/2;
[cartMocapAlignment] = getCartPhasespaceFrameAlignment(alignmentC3D, alignmentRST, ...
  'verbose', 20, 'forceReprocess', forceReprocessEverything, 'phasespaceCalibrationIndeces', phasespaceCalibrationIndeces);

%% Transform data to world frame and generate eye candy
% transform mocap and cart record to world frame using an alignment
[transformedAlignmentPhasespaceRecord] = transformPhasespaceRecordIntoWorldWithMobileCartAlignment(...
  alignmentPhasespaceRecord, alignmentMobileCartRecord, cartMocapAlignment);
[transformedLongWalkPhasespaceRecord] = transformPhasespaceRecordIntoWorldWithMobileCartAlignment(...
  longWalkPhasespaceRecord, longWalkMobileCartRecord, cartMocapAlignment);

% look at mocap markers in world
if (demoMode)
  figure;
  visualizePhasespaceRecord(transformedAlignmentPhasespaceRecord, 'recordName', 'long walk trial');
  % look at cart and mocap markers in world
  figure;
  visualizePhasespaceRecord(transformedAlignmentPhasespaceRecord, 'recordName', 'long walk trial', 'plotDropoutRateHistogram', 0);
  visualizeMobileCartMotion(alignmentMobileCartMotion, 'recordName', 'alignment trial')
  axis equal;
  % make video of mocap markers and cart in world, stationary viewpoint
  figure;
  makeVideoOfPhasespaceRecord('test.avi', transformedAlignmentPhasespaceRecord)
  makeVideoOfPhasespaceRecord('test2.avi', transformedAlignmentPhasespaceRecord, ...
    'mobileCartMotion', alignmentMobileCartMotion);
end

% figure;
% visualizeMobileCartMotion(longWalkMobileCartMotion, 'recordName', 'long walk trial', ...
%   'orientationPlotInterval', 100)
% figure
% makeVideoOfPhasespaceRecord('test.avi', transformedLongWalkPhasespaceRecord, ...
%   'mobileCartMotion', longWalkMobileCartRecord);


% stretchStepsPhasespaceRecord = loadPhasespaceRecord([dataDir 'stretchSteps-01.c3d'], ...
%   'forceReprocess', forceReprocessPhasespace);
% makeVideoOfPhasespaceRecord('test.avi', stretchStepsPhasespaceRecord);




% video, track target, dolly camera

obstaclesPhasespaceRecord = loadPhasespaceRecord([dataDir 'obstacles-01.c3d'], ...
  'forceReprocess', forceReprocessPhasespace);
obstaclesCartFile = [dataDir 'blastobstacles1.rst'];
[transformedObstaclesPhasespaceRecord] = transformPhasespaceRecordIntoWorldWithMobileCartAlignment(...
  obstaclesPhasespaceRecord, obstaclesCartFile, cartMocapAlignment);
makeVideoOfPhasespaceRecord('test.avi', transformedObstaclesPhasespaceRecord, ...
  'mobileCartMotion', obstaclesCartFile, 'trackMarkers', 1);





















