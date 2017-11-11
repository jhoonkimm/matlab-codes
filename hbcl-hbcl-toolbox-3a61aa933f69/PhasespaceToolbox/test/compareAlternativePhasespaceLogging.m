
%%
% customDirectoryName = ...
%   'C:\Users\jrebula\myProjects\IMUGaitAnalysis\rawData\unstructuredTrials\20110509-testCustomCode\phasespaceCustomLog1Fixed.txt';
% directoryName = 'C:\Users\jrebula\myProjects\IMUGaitAnalysis\rawData\unstructuredTrials\20110509-testCustomCode\';

% customDirectoryName = ...
%   'C:\Users\jrebula\myProjects\IMUGaitAnalysis\rawData\unstructuredTrials\20110510-testCustomCode\noFlushing\phasespaceCustomLog-testCustomCode1.txt';
% directoryName = ...
%   'C:\Users\jrebula\myProjects\IMUGaitAnalysis\rawData\unstructuredTrials\20110510-testCustomCode\noFlushing\';

% directoryName = ...
%   'C:\Users\jrebula\myProjects\IMUGaitAnalysis\rawData\unstructuredTrials\20110510-testCustomCode\flushing\';
% customDirectoryName = ...
%   [directoryName '\' 'phasespaceCustomLog1-testCustomCodeWithFlushing1.txt'];


% directoryName = ...
%   'C:\Users\jrebula\myProjects\IMUGaitAnalysis\rawData\unstructuredTrials\20110512-customCodeSyncTest\';
% customDirectoryName = ...
%   [directoryName 'phasespaceCustomLog1-testTreadmillSync.txt'];
% forcePlateFilename = ...
%   [directoryName '20110512-treadmillLogFile.xlsx'];


% directoryName = ...
%   'C:\Users\jrebula\myProjects\IMUGaitAnalysis\rawData\unstructuredTrials\20110512-customCodeSyncTest\test2\';
% customDirectoryName = ...
%   ...   [directoryName 'phasespaceCustomLog3-testTreadmillSync2.txt'];
%   [directoryName 'phasespaceCustomLog2.txt'];
% forcePlateFilename = ...
%   [directoryName 'testPhasespaceSync2.xlsx'];


directoryName = ...
  'C:\Users\jrebula\myProjects\IMUGaitAnalysis\rawData\unstructuredTrials\20110512-customCodeSyncTest\test3\';
customDirectoryName = ...
  [directoryName 'phasespaceCustomLog1.txt'];
forcePlateFilename = ...
  [directoryName 'testPhasespaceSync3-1.xlsx'];

compareC3d = 0;

%%
customPhasespaceLog = processCustomPhasespaceLog(customDirectoryName);
if (compareC3d)
  standardPhasespaceLog = processPhasespaceTrialIntoStructure(directoryName);
end

%%
% markerNames = {'76', ...
%     '77', ...
%     '78', ...
%     '79', ...
%     '80', ...
%     '81', ...
%     '82', ...
%     '83'};
  

% markerNames = {'68', ...
%     '69', ...
%     '70', ...
%     '71'};

markerNames = {'66'}; %'4'}; %

numCols = 1;
for i = 1:length(markerNames)
  markerName = markerNames{i};
  
  %     plot(customPhasespaceLog.(['M' markerName]));
  customMarkerData = customPhasespaceLog.(['M' markerName]);
  customMarkerData(:,1) = -1 * customMarkerData(:,1);
  
  goodIndecesCustom = ~isnan(customMarkerData(:,1));
  
  if (compareC3d)
    standardMarkerData = standardPhasespaceLog.phasespaceData.(['M00' markerName]);
    
    subplot(length(markerNames) / numCols, numCols, i);
    plot(zeros(1, length(standardMarkerData)), 'k--');
    hold on;
    plot(zeros(1, length(standardMarkerData)), 'k');
    
    plot(standardMarkerData, '--', 'LineWidth', 2);
    hold on;
    
    goodIndecesStandard = ~isnan(standardMarkerData(:,1));
    markerDataSplinedStandard = spline(find(goodIndecesStandard), ...
      standardMarkerData(goodIndecesStandard,1), 1:length(standardMarkerData(:,1)));
    
    
    markerDataSplinedCustom = spline(find(goodIndecesCustom), ...
      customMarkerData(goodIndecesCustom,1), 1:length(customMarkerData(:,1)));
    
    %     [cors, lags] = xcorr(customMarkerData(goodIndecesCustom, 1) - mean(customMarkerData(goodIndecesCustom, 1)), ...
    %       standardMarkerData(goodIndecesStandard, 1) - mean(standardMarkerData(goodIndecesStandard, 1)));
    [cors, lags] = xcorr(markerDataSplinedCustom - mean(markerDataSplinedCustom), ...
      markerDataSplinedStandard - mean(markerDataSplinedStandard));
    [temp, bestLagIndex] = max(cors);
    fprintf('time lag from xcorr = %gs or %g phasespace frames\n', lags(bestLagIndex) / 480, lags(bestLagIndex));
  end


  customMarkerOffset = 0; %-lags(bestLagIndex); %(lags(bestLagIndex) + 0 * length(customMarkerData)); %-4000;
  customMarkerIndeces = (1:length(customMarkerData)) + customMarkerOffset;
  
  plot(customMarkerIndeces, customMarkerData);
  hold on;
  
  dataLengthDifference = length(customMarkerData) - length(standardMarkerData);
  
  title(['M' markerName sprintf(' data length difference = %g', dataLengthDifference)]);
end
if (compareC3d)
  legend('c3d file', 'custom log file');
else
  legend('custom log file');
end

%%
% now compare the custom phasespace log with the labview force plate
% recording:
forcePlateData = xlsread(forcePlateFilename, 2);

%%
ledVoltageLevel = forcePlateData(:, 2); %14); %13); %

voltageIsLow = ledVoltageLevel < 3;
markerIsNan = isnan(customMarkerData(:,1));

figure
ax(1) = subplot(2, 1, 1);
voltageTime = (0:(length(ledVoltageLevel) - 1)) * (1/(960*2));
plot(voltageTime, ledVoltageLevel);
hold on;
title('led strobe voltage from labview log');
markerTime = (0:(length(markerIsNan) - 1)) * (1/480);
ax(2) = subplot(2, 1, 2);
plot(markerTime, markerIsNan, 'k', 'LineWidth', 3);
linkaxes(ax,'x');
title('marker is NaN in custom log');

markerIsNanAtVoltageTimes = spline(markerTime, markerIsNan, voltageTime) > 0.99;
voltagesWhenMarkerIsNan = ledVoltageLevel;
voltagesWhenMarkerIsNan(~markerIsNanAtVoltageTimes) = NaN;

subplot(2, 1, 1);
plot(voltageTime, voltagesWhenMarkerIsNan, 'LineWidth', 2, 'color', 'b'); %voltageIsLow);
legend('strobe voltage', 'strobe voltage when marker is NaN');




