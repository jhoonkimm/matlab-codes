
%% found problem... I was assuming 100hz, huge problem, we need to use the
% measured clock ticks from the stellaris for the sample times.

close all;
% clear all;
clc;

forceReprocessEverything = 0; %1; %
dataDir = 'C:\Users\jrebula\myProjects\IMUGaitAnalysis\rawData\overgroundVariabilityExperiment\20110512-PS and Cart Data_larkinsh\calibration\';
alignmentC3D = [dataDir '20110512_calibration-01.c3d'];
alignmentRST = [dataDir 'cartData_20110512_calibration1.rst'];

% dataDir = ...
% 'C:\Users\jrebula\myProjects\IMUGaitAnalysis\rawData\overgroundVariabilityExperiment\20110519-PS and Cart Data_caoke\Calibration\';
% alignmentC3D = [dataDir ''];
% alignmentRST = [dataDir ''];


previousCalibration = load([dataDir 'optimizedCalibrationTransformation.mat']);

% previousCalibrationCombinedFile = load([dataDir 'calibrationOriginal.mat']);
previousCalibrationCombinedFile = load([dataDir 'calibrationEditted.mat']);

leftEncoder = previousCalibrationCombinedFile.phasespaceTrial.cartSensors.leftEncoder;
rightEncoder = previousCalibrationCombinedFile.phasespaceTrial.cartSensors.rightEncoder;

% cartTimes = (0:(length(leftEncoder)-1)) * 0.01;
plot(cartTimes, leftEncoder, cartTimes, rightEncoder);
hold on

marker = previousCalibrationCombinedFile.phasespaceTrial.phasespaceData.M004;
phasespaceTimes = (0:(length(marker)-1)) * 1/480;
plot(phasespaceTimes, marker * 1, 'k');





return;
%%
forceReprocessPhasespace = 0 || forceReprocessEverything;
alignmentPhasespaceRecord = loadPhasespaceRecord(alignmentC3D, 'verbose', 0, 'forceReprocess', forceReprocessPhasespace);
% visualizePhasespaceRecord(alignmentPhasespaceRecord, 'plotInterval', 5, 'recordName', 'alignment trial');

alignmentMobileCartRecord = loadMobileCartRecord(alignmentRST, 'verbose', 0, 'forceReprocess', forceReprocessEverything);
[alignmentMobileCartMotion] = integrateMobileCartMotion(alignmentMobileCartRecord);

% phasespaceCalibrationIndeces = []; %1:length(alignmentPhasespaceRecord.(alignmentPhasespaceRecord.markerNames{1}))/2;
% [cartMocapAlignment] = getCartPhasespaceFrameAlignment(alignmentC3D, alignmentRST, ...
%   'verbose', 20, 'forceReprocess', forceReprocessEverything, 'phasespaceCalibrationIndeces', phasespaceCalibrationIndeces);


%%
metersPerEncoderCount = 2.1494e-04;
cartMocapAlignment = ...
  [previousCalibration.frameAlignment.phasespaceScale / (metersPerEncoderCount * 1000); ...
  previousCalibration.frameAlignment.transformationFromPhasespaceToCart / 1000];

[transformedAlignmentPhasespaceRecord] = transformPhasespaceRecordIntoWorldWithMobileCartAlignment(...
  alignmentPhasespaceRecord, alignmentMobileCartRecord, cartMocapAlignment);

figure;
visualizePhasespaceRecord(transformedAlignmentPhasespaceRecord, 'recordName', 'long walk trial', 'plotDropoutRateHistogram', 0);


return
%%

[transformedAlignmentPhasespaceRecord] = transformPhasespaceRecordIntoWorldWithMobileCartAlignment(...
  alignmentPhasespaceRecord, alignmentMobileCartRecord, cartMocapAlignment);

figure;
visualizePhasespaceRecord(transformedAlignmentPhasespaceRecord, 'recordName', 'long walk trial', 'plotDropoutRateHistogram', 0);
visualizeMobileCartMotion(alignmentMobileCartMotion, 'recordName', 'alignment trial')
axis equal;
figure;
makeVideoOfPhasespaceRecord('test2.avi', transformedAlignmentPhasespaceRecord, ...
  'mobileCartMotion', alignmentMobileCartMotion);

return;



