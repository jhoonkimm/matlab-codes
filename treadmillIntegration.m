
clear all; close all; clc;

% Treadmill Integration for Jay's study
% November 3rd 2017
%---------------------------------------------------------------------------------------------------------------------------------------%

folder = 'SBD_001';

customLog = 'phasespaceCustomLog1_11072017-PWS_01.txt';
treadmillData = '11072017-PWS_01.mat';
c3dFile = 'capture-20171107-PWS_01_CLEANED.c3d';

%---------------------------------------------------------------------------------------------------------------------------------------%
path = 'C:\Users\hbclStudent\Desktop\Jay';
list = {'75Baseline','150Baseline','PWSBaseline','Split','Post'};

customPhasespaceLogName = fullfile(path,folder,customLog);
savedTreadmillDataFilename = fullfile(path,folder,treadmillData);
c3dFilenameToLoad = fullfile(path,folder,c3dFile);

outputFilename = sprintf('integrated_%s_%s',folder,c3dFile);

[ phasespaceAndTreadmillStructure ] = savePhasespaceAndTreadmillDataInC3D(...
  customPhasespaceLogName, savedTreadmillDataFilename, c3dFilenameToLoad, outputFilename,...
  'c3dUnitsAreInMeters',1,...
  'shouldFillGapsInMarkerData',0);