
%%
%
% test on a bit of real data
%

% clear all;
close all;
format compact;

addpath(genpath('C:\\Users\\lafrost\\Desktop\\SvnRepositories\\DataAnalysis\\branches\\HBCL Gait Analysis Toolbox v1.3\\PhasespaceToolbox\\'));

for trialNumber = 2:27
    directory = 'C:\\Users\\lafrost\\Desktop\\08.22.2013\\';
    filename = [directory sprintf('Trial-%02d.c3d', trialNumber)];
    % filledFilename = [directory 'Trial-02 filled.c3d'];
    standingFilename = [directory 'Trial-01.c3d'];
    templateFilename = [directory 'Model_NoTroc.mdh'];
    
    %%
    autofilledFilename = [directory sprintf('Trial-%02d noGTModel autoFilled.c3d', trialNumber)];
    if (~exist(autofilledFilename, 'file'))
        [allInterpolatedMarkers] = autoTrackC3DFile(filename, standingFilename, templateFilename);
        % allInterpolatedMarkers = [];
        % allInterpolatedMarkers.M00 = zeros(2^16 * 4, 3);
        % allInterpolatedMarkers
        writeC3DFile(autofilledFilename, allInterpolatedMarkers);
    end
    
    
end
%%
% trackedC3dTrial = loadPhasespaceRecord(autofilledFilename, ...
%   'c3dUnitsAreInMeters', 1, 'forceReprocess', 0);
%%
% size(allInterpolatedMarkers.M00)
% size(trackedC3dTrial.M00)
% size(allInterpolatedMarkers.M00) - size(trackedC3dTrial.M00)



