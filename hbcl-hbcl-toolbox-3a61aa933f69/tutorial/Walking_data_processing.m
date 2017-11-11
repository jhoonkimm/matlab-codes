clear all
close all
%% 
addpath(genpath('..'))

load('normal.mat')
v=-1.40;
cut_off_freq = 25;
Force_freq = 960;
Marker_freq = 480;
N=1000;
%%
i=1;
 [RHS LHS RTO LTO]=findEventIndicesForWalkingOnTreadmill(RGRF{1},LGRF{1},Force_freq);
 [GoodR, GoodL] = removeBadStridesFromWalkingOnTreadmill(RGRF{1},LGRF{1},RHS,LHS,RTO,LTO,'outputValues','goodStrideNumbers');
    samples = size(RGRF{1},1);
    PLank{i} = interpGaitCycle(PLANK{1}, samples);
    PLkne{i} = interpGaitCycle(PLKNE{1}, samples);
    PLhip{i} = interpGaitCycle(PLHIP{1}, samples);
    PRank{i} = interpGaitCycle(PRANK{1}, samples);
    PRkne{i} = interpGaitCycle(PRKNE{1}, samples);
    PRhip{i} = interpGaitCycle(PRHIP{1}, samples);
    MLank{i} = interpGaitCycle(MLANK{1}, samples);
    MLkne{i} = interpGaitCycle(MLKNE{1}, samples);
    MLhip{i} = interpGaitCycle(MLHIP{1}, samples);
    MRank{i} = interpGaitCycle(MRANK{1}, samples);
    MRkne{i} = interpGaitCycle(MRKNE{1}, samples);
    MRhip{i} = interpGaitCycle(MRHIP{1}, samples);
    ALank{i} = interpGaitCycle(ALANK{1}, samples);
    ALkne{i} = interpGaitCycle(ALKNE{1}, samples);
    ALhip{i} = interpGaitCycle(ALHIP{1}, samples);
    ARank{i} = interpGaitCycle(ARANK{1}, samples);
    ARkne{i} = interpGaitCycle(ARKNE{1}, samples);
    ARhip{i} = interpGaitCycle(ARHIP{1}, samples);
    %%
     for j = 1 : length(LHS)-1
        tL{i}(j) = (LHS(j+1)-LHS(j))/Force_freq;
        RGRFzR{i}(:,j) = interpGaitCycle(RGRF{i}(LHS(j):LHS(j+1),3), N);
        RGRFyR{i}(:,j) =- interpGaitCycle(RGRF{i}(LHS(j):LHS(j+1),2), N);
        RGRFxR{i}(:,j) = interpGaitCycle(RGRF{i}(LHS(j):LHS(j+1),1), N);
        LGRFzR{i}(:,j) = interpGaitCycle(LGRF{i}(LHS(j):LHS(j+1),3), N);
        LGRFyR{i}(:,j) = -interpGaitCycle(LGRF{i}(LHS(j):LHS(j+1),2), N);
        LGRFxR{i}(:,j) = interpGaitCycle(LGRF{i}(LHS(j):LHS(j+1),1), N);
        
        PRank_gait{i}(:,j) = interpGaitCycle(PRank{i}(LHS(j):LHS(j+1),1), N);
        PRkne_gait{i}(:,j) = interpGaitCycle(PRkne{i}(LHS(j):LHS(j+1),1), N);
        PRhip_gait{i}(:,j) = interpGaitCycle(PRhip{i}(LHS(j):LHS(j+1),1), N);
        PLank_gait{i}(:,j) = interpGaitCycle(PLank{i}(LHS(j):LHS(j+1),1), N);
        PLkne_gait{i}(:,j) = interpGaitCycle(PLkne{i}(LHS(j):LHS(j+1),1), N);
        PLhip_gait{i}(:,j) = interpGaitCycle(PLhip{i}(LHS(j):LHS(j+1),1), N);
        
        MRank_gait{i}(:,:,j) = interpGaitCycle(MRank{i}(LHS(j):LHS(j+1),:), N);
        MRkne_gait{i}(:,:,j) = interpGaitCycle(MRkne{i}(LHS(j):LHS(j+1),:), N);
        MRhip_gait{i}(:,:,j) = interpGaitCycle(MRhip{i}(LHS(j):LHS(j+1),:), N);
        MLank_gait{i}(:,:,j) = interpGaitCycle(MLank{i}(LHS(j):LHS(j+1),:), N);
        MLkne_gait{i}(:,:,j) = interpGaitCycle(MLkne{i}(LHS(j):LHS(j+1),:), N);
        MLhip_gait{i}(:,:,j) = interpGaitCycle(MLhip{i}(LHS(j):LHS(j+1),:), N);
        
        ARank_gait{i}(:,:,j) = interpGaitCycle(ARank{i}(LHS(j):LHS(j+1),:), N);
        ARkne_gait{i}(:,:,j) = interpGaitCycle(ARkne{i}(LHS(j):LHS(j+1),:), N);
        ARhip_gait{i}(:,:,j) = interpGaitCycle(ARhip{i}(LHS(j):LHS(j+1),:), N);
        ALank_gait{i}(:,:,j) = interpGaitCycle(ALank{i}(LHS(j):LHS(j+1),:), N);
        ALkne_gait{i}(:,:,j) = interpGaitCycle(ALkne{i}(LHS(j):LHS(j+1),:), N);
        ALhip_gait{i}(:,:,j) = interpGaitCycle(ALhip{i}(LHS(j):LHS(j+1),:), N);
      
        [Lank_work{i}(j) Lank_work_neg{i}(j)] = findPosNegWork(PLank_gait{i}(:,j),tL{i}(j));
        [Lkne_work{i}(j) Lkne_work_neg{i}(j)] = findPosNegWork(PLkne_gait{i}(:,j),tL{i}(j));
        [Lhip_work{i}(j) Lhip_work_neg{i}(j)]= findPosNegWork(PLhip_gait{i}(:,j),tL{i}(j));
        
        [Rank_work{i}(j) Rank_work_neg{i}(j)] = findPosNegWork(PRank_gait{i}(:,j),tL{i}(j));
        [Rkne_work{i}(j) Rkne_work_neg{i}(j)] = findPosNegWork(PRkne_gait{i}(:,j),tL{i}(j));
        [Rhip_work{i}(j) Rhip_work_neg{i}(j)] = findPosNegWork(PRhip_gait{i}(:,j),tL{i}(j));
        
        Lank_workrate{i}(j) = Lank_work{i}(j)*2/tL{i}(j);
        Lkne_workrate{i}(j) = Lkne_work{i}(j)*2/tL{i}(j);
        Lhip_workrate{i}(j) = Lhip_work{i}(j)*2/tL{i}(j);
        Rank_workrate{i}(j) = Rank_work{i}(j)*2/tL{i}(j);
        Rkne_workrate{i}(j) = Rkne_work{i}(j)*2/tL{i}(j);
        Rhip_workrate{i}(j) = Rhip_work{i}(j)*2/tL{i}(j);
       

     end
       step_frequency_ave(i)=mean(2./tL{i}(GoodL));
    step_frequency_std(i)=std(2./tL{i}(GoodL));
    RGRF_aveR(:,1,i)=mean(RGRFxR{i}(:,GoodL),2);
    RGRF_aveR(:,2,i)=-mean(RGRFyR{i}(:,GoodL),2);
    RGRF_aveR(:,3,i)=mean(RGRFzR{i}(:,GoodL),2);
    LGRF_aveR(:,1,i)=mean(LGRFxR{i}(:,GoodL),2);
    LGRF_aveR(:,2,i)=-mean(LGRFyR{i}(:,GoodL),2);
    LGRF_aveR(:,3,i)=mean(LGRFzR{i}(:,GoodL),2);
    
    PRank_gait_ave(:,i) = mean(PRank_gait{i}(:,GoodL),2);
    PRkne_gait_ave(:,i) = mean(PRkne_gait{i}(:,GoodL),2);
    PRhip_gait_ave(:,i) = mean(PRhip_gait{i}(:,GoodL),2);
    PLank_gait_ave(:,i) = mean(PLank_gait{i}(:,GoodL),2);
    PLkne_gait_ave(:,i) = mean(PLkne_gait{i}(:,GoodL),2);
    PLhip_gait_ave(:,i) = mean(PLhip_gait{i}(:,GoodL),2);

    MRank_gait_ave(:,:,i) = mean(MRank_gait{i}(:,:,GoodL),3);
    MRkne_gait_ave(:,:,i) = mean(MRkne_gait{i}(:,:,GoodL),3);
    MRhip_gait_ave(:,:,i) = mean(MRhip_gait{i}(:,:,GoodL),3);
    MLank_gait_ave(:,:,i) = mean(MLank_gait{i}(:,:,GoodL),3);
    MLkne_gait_ave(:,:,i) = mean(MLkne_gait{i}(:,:,GoodL),3);
    MLhip_gait_ave(:,:,i) = mean(MLhip_gait{i}(:,:,GoodL),3);

    ARank_gait_ave(:,:,i) = mean(ARank_gait{i}(:,:,GoodL),3);
    ARkne_gait_ave(:,:,i) = mean(ARkne_gait{i}(:,:,GoodL),3);
    ARhip_gait_ave(:,:,i) = mean(ARhip_gait{i}(:,:,GoodL),3);
    ALank_gait_ave(:,:,i) = mean(ALank_gait{i}(:,:,GoodL),3);
    ALkne_gait_ave(:,:,i) = mean(ALkne_gait{i}(:,:,GoodL),3);
    ALhip_gait_ave(:,:,i) = mean(ALhip_gait{i}(:,:,GoodL),3);

   [RComPR{i} LComPR{i} RComP_aveR{i} LComP_aveR{i} RComP_stdR{i} LComP_stdR{i} velocityR{i}]=...
        CalComPower (RGRFxR{i}, RGRFyR{i}, RGRFzR{i},LGRFxR{i},LGRFyR{i},LGRFzR{i},tL{i}',GoodL,v);
for j = 1 : size(RComPR{i},2) 
        [HS{i}(j) RB{i}(j) PL{i}(j) PO{i}(j)]=CalPhaseWork(LComPR{i}(:,j),LComPR{i}(:,j),tL{i}(GoodL(j)));
       
end

%%
figure(1);clf;hold all;
gait_percent = 0.1:0.1:100;
    RCOM = RComP_aveR{1} ;
    LCOM = LComP_aveR{1} ;
    
% plot(gait_percent,RCOM,'linewidth', 2);
plot(gait_percent,LCOM,'-','linewidth', 2);
hfill=fill([gait_percent gait_percent(end:-1:1)]',[LComP_aveR{1}+LComP_stdR{1} ;LComP_aveR{1}(end:-1:1)-2*LComP_stdR{1}(end:-1:1)],'r');
set(hfill,'FaceAlpha',0.5,'EdgeColor','r')

line([0 100],[0 0],'color','k');
xlabel('% gait cycle')
ylabel('COM work rate')

%%

AnkleAngles = -ALank_gait_ave(:,1)+65;
AnkleMoment = MLank_gait_ave(:,1);
AnklePower = PLank_gait_ave;

KneeAngles = ALkne_gait_ave(:,1);
KneeMoment = MLkne_gait_ave(:,1);
KneePower = PLkne_gait_ave;

HipAngles = ALhip_gait_ave(:,1);
HipMoment = MLhip_gait_ave(:,1);
HipPower = PLhip_gait_ave;
plot3x3(AnkleAngles,AnkleMoment,AnklePower,KneeAngles,KneeMoment,KneePower,HipAngles,HipMoment,HipPower,'powerUnits','W','momentUnits','N-m');

%%


