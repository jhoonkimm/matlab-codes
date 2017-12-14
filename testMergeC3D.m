
% customLog = 'phasespaceCustomLog6_11142017_Split_01.txt';
% treadmillData = '11142017_Split_01';

c3dFile1 = 'capture-20171121-Split-01_001-CLEANED.c3d';
c3dFile2 = 'capture-20171121-Split-01_002-CLEANED.c3d';
c3dFile3 = 'capture-20171121-Split-01_003-CLEANED.c3d';

c3dFileCombined = 'capture-20171121-Split-01_combined-CLEANED.c3d';

%---------------------------------------------------------------------------------------------------------------------------------------%
path = 'E:\RESEARCH\MDSC 508';
folder = 'SBD_004';
% list = {'75Baseline','150Baseline','PWSBaseline','Split','Post'};

% customPhasespaceLogName = fullfile(path,folder,customLog);
% savedTreadmillDataFilename = fullfile(path,folder,treadmillData);

c3dFilename1 = fullfile(path,folder,c3dFile1);
c3dFilename2 = fullfile(path,folder,c3dFile2);
c3dFilename3 = fullfile(path,folder,c3dFile3);
c3dFilenameToSave = fullfile(path,folder,c3dFileCombined);

% outputFilename = sprintf('integrated_%s_%s',folder,c3dFile);

mergeC3D(c3dFilename1, c3dFilename2, c3dFilename3, c3dFilenameToSave)

%%
% 
% [ phasespaceAndTreadmillStructure ] = savePhasespaceAndTreadmillDataInC3D2018(...
%       fullfile(path,folder,'phasespaceCustomLog6_11212017_Split_01.txt'), fullfile(path,folder,'11212017_Split_01.mat'), fullfile(path,folder,'capture-20171121-Split-01_combined-CLEANED'), fullfile(path,folder,sprintf('integrated_%s','capture-20171121-Split-01_combined-CLEANED')),...
%       'c3dUnitsAreInMeters',1,...
%       'shouldFillGapsInMarkerData',0)
  
[ phasespaceAndTreadmillStructure ] = savePhasespaceAndTreadmillDataInC3D2018(...
      fullfile(path,folder,'phasespaceCustomLog3_11212017_PWSBaseline_02.txt'), ...
      fullfile(path,folder,'11212017_PWSBaseline_02.mat'),...
      fullfile(path,folder,'capture-20171121-PWSBaseline-02-CLEANED.c3d'),...
      fullfile(path,folder,sprintf('integrated_%s','capture-20171121-PWSBaseline-02-CLEANED')),...
      'c3dUnitsAreInMeters',1,...
      'shouldFillGapsInMarkerData',0)
  
  %%
  load('E:\RESEARCH\MDSC 508\SBD_004\sbd004powers.mat')
  l_ank_power = l_ank_power{1};
  r_ank_power = r_ank_power{1};
  l_kne_power = l_kne_power{1};
  r_kne_power = r_kne_power{1};
  l_hip_power = l_hip_power{1};
  r_hip_power = r_hip_power{1};
  al_ank_power = mean(l_ank_power,2);
  ar_ank_power = mean(r_ank_power,2);
  al_kne_power = mean(l_kne_power,2);
  ar_kne_power = mean(r_kne_power,2);
  al_hip_power = mean(l_hip_power,2);
  ar_hip_power = mean(r_hip_power,2);
    %%
  load('E:\RESEARCH\MDSC 508\SBD_004\sbd004powerspws.mat')
  l_ank_powerpws = l_ank_power{1};
  r_ank_powerpws = r_ank_power{1};
  l_kne_powerpws = l_kne_power{1};
  r_kne_powerpws = r_kne_power{1};
  l_hip_powerpws = l_hip_power{1};
  r_hip_powerpws = r_hip_power{1};
  al_ank_powerpws = mean(l_ank_powerpws,2);
  ar_ank_powerpws = mean(r_ank_powerpws,2);
  al_kne_powerpws = mean(l_kne_powerpws,2);
  ar_kne_powerpws = mean(r_kne_powerpws,2);
  al_hip_powerpws = mean(l_hip_powerpws,2);
  ar_hip_powerpws = mean(r_hip_powerpws,2);
  
 %%
 
FM = loadForcesFromHBCLBertecTreadmillMatFile('E:\RESEARCH\MDSC 508\SBD_004\11212017_Split_01.mat');


LGRF = FM.left.groundReactionForces;
RGRF = FM.right.groundReactionForces;
LGRM = FM.left.groundReactionMoments;
RGRM = FM.right.groundReactionMoments;
[RHS,LHS,RTO,LTO] = findEventIndicesForWalkingOnTreadmill(RGRF,LGRF);
%%
FMpws = loadForcesFromHBCLBertecTreadmillMatFile('E:\RESEARCH\MDSC 508\SBD_004\11212017_PWSBaseline_02.mat');
LGRFpws = FMpws.left.groundReactionForces;
RGRFpws = FMpws.right.groundReactionForces;
LGRMpws = FMpws.left.groundReactionMoments;
RGRMpws = FMpws.right.groundReactionMoments;
[RHSpws,LHSpws,RTOpws,LTOpws] = findEventIndicesForWalkingOnTreadmill(RGRFpws,LGRFpws);


%%
orgR = []; orgL = [];

for k = 2:length(RHS)
    idx = unique(round((RHS(k-1):RHS(k))/2));
    orgR = [orgR,interpGaitCycle(ar_kne_power(idx),500)];
    orgL = [orgL,interpGaitCycle(al_kne_power(idx),500)];
end
%%
orgRpws = []; orgLpws = [];
hold on
axis([1 500 -0.1 0.1])
for k = 2:length(RHSpws)
    idx = unique(round((RHSpws(k-1):RHSpws(k))/2));
    orgRpws = [orgRpws,interpGaitCycle(ar_ank_powerpws(idx),500)];
    orgLpws = [orgLpws,interpGaitCycle(al_kne_powerpws(idx),500)];
    plot(1:500,orgRpws(:,k-1))
pause(0.05)
end

hold off
%%

plot(1:500,mean(orgR,2),'b-',1:500,mean(orgL,2),'r-',1:500,mean(orgRpws,2),'b--',1:500,mean(orgLpws,2),'r--')
  