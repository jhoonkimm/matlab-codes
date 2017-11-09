clc;

metFilename = {'C:\Users\hbclStudent\Desktop\Jay\MATLAB Scripts\SBD_001\10146402_resting.txt'
               'C:\Users\hbclStudent\Desktop\Jay\MATLAB Scripts\SBD_001\10146402_pws.txt'
               'C:\Users\hbclStudent\Desktop\Jay\MATLAB Scripts\SBD_001\10146402_150baseline.txt'
               'C:\Users\hbclStudent\Desktop\Jay\MATLAB Scripts\SBD_001\10146402_75baseline.txt'
               'C:\Users\hbclStudent\Desktop\Jay\MATLAB Scripts\SBD_001\10146402_split.txt'};
tot = [];
% for i = 1:length(metFilename)
    metOutput1 = oxyconProcessing(metFilename{i},'brockway',1);
    metOutput2 = oxyconProcessing(metFilename{i},'brockway',0,'massSpecific',0);
%     tot = [tot;metOutput];
% end

%%
clear all; close all; clc;

matFilename = 'E:\RESEARCH\MDSC 528\TEST_SBD_001\logfile_75Baseline_01.mat';
forceOutput  = loadForcesFromHBCLBertecTreadmillMatFile(matFilename,'shouldFilter',1);

%groundReactionForces
% 1st col = x or lateral
% 2nd col = y or fore-aft
% 3rd col = z or vertical

RGRF = forceOutput.right.groundReactionForces;
LGRF = forceOutput.left.groundReactionForces;
RGRM = forceOutput.right.groundReactionMoments;
LGRM = forceOutput.left.groundReactionMoments;

[RHS,LHS,RTO,LTO] = findEventIndicesForWalkingOnTreadmill(RGRF,LGRF, 'dataFrequency', 960,'showPlots',false);

%% Finding COP per cycle

%each Rcyc/Lcyc has
% columns
%1 = vert force, 2 = x GRM, 3 = y GRM, 4 = x COP, 4 = y COP

hold on
for i = 1:length(RHS)
    Rcyc{i,1} = interpGaitCycle(RGRF(RHS(i):RTO(i),3));
    Lcyc{i,1} = interpGaitCycle(LGRF(LHS(i):LTO(i),3));
    Rcyc{i,2} = interpGaitCycle(RGRM(RHS(i):RTO(i),1));
    Lcyc{i,2} = interpGaitCycle(LGRM(LHS(i):LTO(i),1));
    Rcyc{i,3} = interpGaitCycle(RGRM(RHS(i):RTO(i),2));
    Lcyc{i,3} = interpGaitCycle(LGRM(LHS(i):LTO(i),2));
    for j = 1:1000
        Rcyc{i,4}(j,1) = -Rcyc{i,3}(j)/Rcyc{i,1}(j);
        Lcyc{i,4}(j,1) = -Lcyc{i,3}(j)/Lcyc{i,1}(j);
        Rcyc{i,5}(j,1) = Rcyc{i,2}(j)/Rcyc{i,1}(j);
        Lcyc{i,5}(j,1) = Lcyc{i,2}(j)/Lcyc{i,1}(j);
    end
    
    %plot(Rcyc{i,5},Rcyc{i,5},'b-',Lcyc{i,4},Lcyc{i,5},'r-')
    pause(0.05)
%     axis([-0.7 -0.1 0.5 2])
    
end
hold off




