%% plot hip joint power
close all
linewidth=2;

figure(1);
subplot(3,3,7);hold on
plot([0.1:0.1:100]',PRhip_aveR,'linewidth',linewidth);
% axis([0 100 -0.1 0.17]);
% plot([0.1:0.1:100]',PRhips_normal,'linewidth',linewidth,'color',[190 190 190 ]./255);
xlabel('Time(% per stride)','fontsize',10);
ylabel('Power(W/mg(gL)^{0.5})','fontsize',10);
title('Right hip joint power','fontsize',10);
%saveas(11,'Right hip joint power','ai');
%saveas(11,'Right hip joint power.eps','psc2');
%saveas(11,'Right hip joint power','fig');

figure(2);
subplot(3,3,7);hold on
plot([0.1:0.1:100]',PLhip_aveR,'linewidth',linewidth);
% axis([0 100 -0.1 0.17]);
% plot([0.1:0.1:100]',PLhips_normal,'LineWidth',linewidth,'color',[190 190 190 ]./255)
xlabel('Time(% per stride)');
ylabel('Power(W/mg(gL)^{0.5})');
title('Left hip joint power');
% saveas(12,'Left hip joint power','ai');
% saveas(12,'Left hip joint power.eps','psc2');
% saveas(12,'Left hip joint power','fig');

%% Plot knee joint power
figure(1);
subplot(3,3,8);hold on
plot([0.1:0.1:100]',PRkne_aveR,'linewidth',linewidth);
% axis([0 100 -0.1 0.17]);
% plot([0.1:0.1:100]',PRknes_normal,'linewidth',linewidth,'color',[190 190 190 ]./255);
xlabel('Time(% per stride)','fontsize',10);
%ylabel('Power(W/mg(gL)^{0.5})');
title('Right knee joint power','fontsize',10);
%saveas(13,'Right knee joint power','ai');
%saveas(13,'Right knee joint power.eps','psc2');
%saveas(13,'Right knee joint power','fig');

figure(2);
subplot(3,3,8);hold on
plot([0.1:0.1:100]',PLkne_aveR,'linewidth',linewidth);
% axis([0 100 -0.1 0.17]);
% plot([0.1:0.1:100]',PLknes_normal,'linewidth',linewidth,'color',[190 190 190 ]./255);
xlabel('Time(% per stride)');
ylabel('Power(W/mg(gL)^{0.5})');
title('Left knee joint power');
% saveas(14,'Left knee joint power','ai');
% saveas(14,'Left knee joint power','fig');
% saveas(14,'Left knee joint power.eps','psc2');

%% Plot ankle joint power
figure(1);
subplot(3,3,9);hold on
plot([0.1:0.1:100]',PRank_aveR,'linewidth',linewidth);
% plot([0.1:0.1:100]',PRanks_normal,'linewidth',linewidth,'color',[190 190 190 ]./255);
% axis([0 100 -0.1 0.17]);
xlabel('Time(% per stride)','fontsize',10);
%ylabel('Power(W/mg(gL)^{0.5})');
%legend('6.60 kN/m','18.68 kN/m','38.66 kN/m','58.29 kN/m');
title('Right ankle joint power','fontsize',10);
%saveas(15,'Right ankle joint power','ai');
%saveas(15,'Right ankle joint power','fig');
%saveas(15,'Right ankle joint power.eps','psc2');

figure(2);
subplot(3,3,9);hold on
plot([0.1:0.1:100]',PLank_aveR,'linewidth',linewidth);
% plot([0.1:0.1:100]',PLanks_normal,'linewidth',linewidth,'color',[190 190 190 ]./255);
% axis([0 100 -0.1 0.17]);

xlabel('Time(% per stride)');
ylabel('Power(W/mg(gL)^{0.5})');
%legend('6.60 kN/m','18.68 kN/m','38.66 kN/m','58.29 kN/m');
title('Left ankle joint power');
% saveas(16,'Left ankle joint power','ai');
% saveas(16,'Left ankle joint power','fig');
% saveas(16,'Left ankle joint power.eps','psc2');


%% plot hip moment
figure(1);subplot(3,3,4);hold on
% plot([0.1:0.1:100]',MRhips_normal,'linewidth',linewidth,'color',[190 190 190 ]./255);
plot([0.1:0.1:100]',[MRhip_aveR(:,1,1) MRhip_aveR(:,1,2) MRhip_aveR(:,1,3) MRhip_aveR(:,1,4)],'linewidth',linewidth);
% axis([0 100 -0.25 0.1]);
%xlabel('Time(% per stride)');
ylabel('Moment(Nm/mgL)','fontsize',10);
%legend('6.60 kN/m','18.68 kN/m','38.66 kN/m','58.29 kN/m');
title('Right hip joint moment','fontsize',10);
%saveas(17,'Right hip joint moment','ai');
%saveas(17,'Right hip joint moment.eps','psc2');
%saveas(17,'Right hip joint moment','fig');

figure(2);
subplot(3,3,4);hold on
% plot([0.1:0.1:100]',MLhips_normal,'linewidth',linewidth,'color',[190 190 190 ]./255);
plot([0.1:0.1:100]',[MLhip_aveR(:,1,1) MLhip_aveR(:,1,2) MLhip_aveR(:,1,3) MLhip_aveR(:,1,4)],'linewidth',linewidth);
% axis([0 100 -0.25 0.1]);

xlabel('Time(% per stride)');
ylabel('Moment(Nm/mgL)');
%legend('6.60 kN/m','18.68 kN/m','38.66 kN/m','58.29 kN/m');
title('Left hip joint moment');
% saveas(18,'Left hip joint moment','ai');
% saveas(18,'Left hip joint moment','fig');
% saveas(18,'Left hip joint moment.eps', 'psc2')
%% plot knee moment

figure(1);subplot(3,3,5);hold on
plot([0.1:0.1:100]',-[MRkne_aveR(:,1,1) MRkne_aveR(:,1,2) MRkne_aveR(:,1,3) MRkne_aveR(:,1,4)],'linewidth',linewidth);
% plot([0.1:0.1:100]',MRknes_normal,'linewidth',linewidth,'color',[190 190 190 ]./255);
% axis([0 100 -0.25 0.1]);

%xlabel('Time(% per stride)');
%ylabel('Moment(Nm/mgL)');
%legend('6.60 kN/m','18.68 kN/m','38.66 kN/m','58.29 kN/m');
title('Right knee joint moment','fontsize',10);
%saveas(19,'Right knee joint moment','ai');
%saveas(19,'Right knee joint moment','fig');
%saveas(19,'Right knee joint moment.eps', 'psc2')
figure(2);
subplot(3,3,5);hold on
% plot([0.1:0.1:100]',MLknes_normal,'linewidth',linewidth,'color',[190 190 190 ]./255);
plot([0.1:0.1:100]',-[MLkne_aveR(:,1,1) MLkne_aveR(:,1,2) MLkne_aveR(:,1,3) MLkne_aveR(:,1,4)],'linewidth',linewidth);
% axis([0 100 -0.25 0.1]);

xlabel('Time(% per stride)');
ylabel('Moment(Nm//mgL)');
%legend('6.60 kN/m','18.68 kN/m','38.66 kN/m','58.29 kN/m');
title('Left knee joint moment');
% saveas(20,'Left knee joint moment','ai');
% saveas(20,'Left knee joint moment','fig');
% saveas(20,'Left knee joint moment.eps', 'psc2')
clear RRR LLL
%% plot ankle moment

figure(1);subplot(3,3,6);hold on
% plot([0.1:0.1:100]',MRanks_normal,'linewidth',linewidth,'color',[190 190 190 ]./255);

plot([0.1:0.1:100]',[MRank_aveR(:,1,1) MRank_aveR(:,1,2) MRank_aveR(:,1,3) MRank_aveR(:,1,4)],'linewidth',linewidth);
% axis([0 100 -0.25 0.1]);

%xlabel('Time(% per stride)');
%ylabel('Moment(Nm/mgL)');
%legend('6.60 kN/m','18.68 kN/m','38.66 kN/m','58.29 kN/m');
title('Right ankle joint moment','fontsize',10);
name='Right ankle joint moment';
%saveas(21,name,'ai');
%saveas(21,name,'fig');
%saveas(21,[name '.eps'], 'psc2')
figure(2);
subplot(3,3,6);hold on
% plot([0.1:0.1:100]',MLanks_normal,'linewidth',linewidth,'color',[190 190 190 ]./255);
plot([0.1:0.1:100]',[MLank_aveR(:,1,1) MLank_aveR(:,1,2) MLank_aveR(:,1,3) MLank_aveR(:,1,4)],'linewidth',linewidth);
% axis([0 100 -0.25 0.1]);

xlabel('Time(% per stride)','fontsize',10);
ylabel('Moment(Nm/mgL)','fontsize',10);
%legend('6.60 kN/m','18.68 kN/m','38.66 kN/m','58.29 kN/m');
title('Left ankle joint moment','fontsize',10);
name='Left ankle joint moment';
% saveas(22,name,'ai');
% saveas(22,name,'fig');
% saveas(22,[name '.eps'], 'psc2')


%% plot hip ankle

figure(1);subplot(3,3,1);hold on
% plot([0.1:0.1:100]',ARhips_normal,'linewidth',linewidth,'color',[190 190 190 ]./255);
plot([0.1:0.1:100]',[ARhip_aveR(:,1,1) ARhip_aveR(:,1,2) ARhip_aveR(:,1,3) ARhip_aveR(:,1,4) ARhip_aveR(:,1,5)],'linewidth',linewidth);
% axis([0 100 -90 100]);
%xlabel('Time(% per stride)');
ylabel('Degree','fontsize',10);
%legend('6.60 kN/m','18.68 kN/m','38.66 kN/m','58.29 kN/m');
title('Right hip joint angle','fontsize',10);
name='Right hip joint angle';
%saveas(23,name,'ai');
%saveas(23,name,'fig');
%saveas(23,[name '.eps'], 'psc2')
figure(2);
subplot(3,3,1);hold on
% plot([0.1:0.1:100]',ALhips_normal,'linewidth',linewidth,'color',[190 190 190 ]./255);
plot([0.1:0.1:100]',[ALank_aveR(:,1,1) ALank_aveR(:,1,2) ALank_aveR(:,1,3) ALank_aveR(:,1,4)],'linewidth',linewidth);
% axis([0 100 -90 100]);

xlabel('Time(% per stride)','fontsize',10);
ylabel('Degree','fontsize',10);
%legend('6.60 kN/m','18.68 kN/m','38.66 kN/m','58.29 kN/m');
title('Left hip joint angle','fontsize',10);
name='Left hip joint angle';
% saveas(24,name,'ai');
% saveas(24,name,'fig');
% saveas(24,[name '.eps'], 'psc2')

clear RRR LLL
%% plot knee angle

figure(1);subplot(3,3,2);hold on
% plot([0.1:0.1:100]',ARknes_normal,'linewidth',linewidth,'color',[190 190 190 ]./255);
plot([0.1:0.1:100]',-[ARkne_aveR(:,1,1) ARkne_aveR(:,1,2) ARkne_aveR(:,1,3) ARkne_aveR(:,1,4)],'linewidth',linewidth);
axis([0 100 -90 100]);

%xlabel('Time(% per stride)');
%ylabel('Degree');
%legend('6.60 kN/m','18.68 kN/m','38.66 kN/m','58.29 kN/m');
title('Right knee joint angle','fontsize',10);
name='Right knee joint angle';
%saveas(25,name,'ai');
%saveas(25,name,'fig');
%saveas(25,[name '.eps'], 'psc2')
figure(2);
subplot(3,3,2);hold on
% plot([0.1:0.1:100]',ALknes_normal,'linewidth',linewidth,'color',[190 190 190 ]./255);
plot([0.1:0.1:100]',-[ALkne_aveR(:,1,1) ALkne_aveR(:,1,2) ALkne_aveR(:,1,3) ALkne_aveR(:,1,4)],'linewidth',linewidth);
% axis([0 100 -90 100]);

xlabel('Time(% per stride)','fontsize',10);
ylabel('Degree','fontsize',10);
%legend('6.60 kN/m','18.68 kN/m','38.66 kN/m','58.29 kN/m');
title('Left knee joint angle','fontsize',10);
name='Left knee joint angle';
% saveas(26,name,'ai');
% saveas(26,name,'fig');
% saveas(26,[name '.eps'], 'psc2')
clear RRR LLL
%% plot ankle angle

figure(1);subplot(3,3,3);hold on
% plot([0.1:0.1:100]',ARanks_normal,'linewidth',linewidth,'color',[190 190 190 ]./255);
plot([0.1:0.1:100]',[ARank_aveR(:,1,1) ARank_aveR(:,1,2) ARank_aveR(:,1,3) ARank_aveR(:,1,4)],'linewidth',linewidth);
% axis([0 100 -90 100]);

%xlabel('Time(% per stride)');
%ylabel('Degree');
%legend('6.60 kN/m','18.68 kN/m','38.66 kN/m','58.29 kN/m');
title('Right ankle joint angle','fontsize',10);
legend('6.60 kN/m','18.68 kN/m','38.66 kN/m','58.29 kN/m');

name='Right ankle joint angle';
%saveas(27,name,'ai');
%saveas(27,[name '.eps'], 'psc2');
%saveas(27,name,'fig');

figure(2);
subplot(3,3,3);hold on
% plot([0.1:0.1:100]',ALanks_normal,'linewidth',linewidth,'color',[190 190 190 ]./255);
plot([0.1:0.1:100]',[ALank_aveR(:,1,1) ALank_aveR(:,1,2) ALank_aveR(:,1,3) ALank_aveR(:,1,4)],'linewidth',linewidth);
% axis([0 100 -90 100]);

xlabel('Time(% per stride)','fontsize',10);
ylabel('Degree','fontsize',10);
legend('6.60 kN/m','18.68 kN/m','38.66 kN/m','58.29 kN/m');
title('Left ankle joint angle');
name='Left ankle joint angle';

% saveas(28,name,'ai');
% saveas(28,[name '.eps'], 'psc2')
% saveas(28,name,'fig');

saveas(1,'Right3by3','ai');
saveas(1,'Right3by3.eps', 'psc2')
saveas(1,'Right3by3','fig');

saveas(1,'Right3by3','jpg');
saveas(2,'Left3by3','ai');
saveas(2,'Left3by3.eps', 'psc2')
saveas(2,'Left3by3','fig');

saveas(2,'Left3by3','jpg');
%%
figure(3);
plot([0.1:0.1:100]',RFT_ave,'linewidth',linewidth );
xlabel('Time(% per stride)','fontsize',10);
ylabel('Power(W)','fontsize',10);
legend('6.60 kN/m','18.68 kN/m','38.66 kN/m','58.29 kN/m');
title('Right foot inter-segment Power');
name='Right_inter';
saveas(3,name,'ai');
saveas(3,[name '.eps'], 'psc2')
saveas(3,name,'fig');
saveas(3,name,'jpg');
figure(4);
plot([0.1:0.1:100]',LFT_ave,'linewidth',linewidth );
xlabel('Time(% per stride)','fontsize',10);
ylabel('Power(W)','fontsize',10);
legend('6.60 kN/m','18.68 kN/m','38.66 kN/m','58.29 kN/m');
title('Left foot inter-segment Power');
name='Left_inter';
saveas(4,name,'ai');
saveas(4,[name '.eps'], 'psc2')
saveas(4,name,'fig');
saveas(4,name,'jpg');
%%
figure(5);
plot([0.1:0.1:100]',[RComP_aveR{1} RComP_aveR{2} RComP_aveR{3} RComP_aveR{4} RComP_aveR{5}], 'linewidth',linewidth);
hold on;
% plot([0.1:0.1:100]',[LComP_aveR{1} LComP_aveR{2} LComP_aveR{3} LComP_aveR{4}],'--', 'linewidth',linewidth);

xlabel('Time(% per stride)','fontsize',10);
ylabel('Work rate(W)','fontsize',10);
legend('6.60 kN/m','18.68 kN/m','38.66 kN/m','58.29 kN/m');
name='COM work rate';
saveas(5,name,'ai');
saveas(5,[name '.eps'], 'psc2')
saveas(5,name,'fig');
saveas(5,name,'jpg');

figure(6);
plot([0.1:0.1:100]',RGRFzR_ave,'linewidth',linewidth);
xlabel('Time(% per stride)','fontsize',10);
ylabel('Vertical ground reaction force(N)','fontsize',10);
legend('6.60 kN/m','18.68 kN/m','38.66 kN/m','58.29 kN/m');
name='GRFRz';
saveas(6,name,'ai');
saveas(6,[name '.eps'], 'psc2')
saveas(6,name,'fig');
saveas(6,name,'jpg');
