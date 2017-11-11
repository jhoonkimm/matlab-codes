close all
stiffness=[6.6 18.68 38.66 58.29];
%% Plot Collision
figure(1);
errorbar(stiffness,-HS_ave,HS_std);
xlabel('Stiffness(kN/m)');
ylabel('Work(J)');
title('Collision work');
saveas(1,'Collision work','ai');
saveas(1,'Collision work','fig');
saveas(1,'Collision work.eps','psc2');


%% Plot Rebound
figure(2);
errorbar(stiffness,RB_ave,RB_std);
xlabel('Stiffness(kN/m)');
ylabel('Work(J)');
title('Rebound work');
saveas(2,'Rebound work','ai');
saveas(2,'Rebound work','fig');
saveas(2,'Rebound work.eps','psc2');

%% Plot Preload
figure(3);
errorbar(stiffness,-PL_ave,PL_std);
xlabel('Stiffness(kN/m)');
ylabel('Work(J)');
title('Preload work');
saveas(3,'Preload work','ai');
saveas(3,'Preload work','fig');
saveas(3,'Preload work.eps','psc2');

%% Plot Push-off
figure(4);
errorbar(stiffness,PO_ave,PO_std);
xlabel('Stiffness(kN/m)');
ylabel('Work(J)');
title('Push-off work');
saveas(4,'Pushoff work','ai');
saveas(4,'Pushoff work','fig');
saveas(4,'Pushoff work.eps','psc2');

%% Plot Hip work
figure(5);hold on;
errorbar(stiffness,WLhip_pos_ave,WLhip_pos_std);
xlabel('Stiffness(kN/m)');
ylabel('Work(J)');
title('Left hip joint work');
errorbar(stiffness,WLhip_neg_ave,WLhip_neg_std,'r');
hold off
saveas(5,'Left hip joint work','ai');
saveas(5,'Left hip joint work.eps','psc2');
saveas(5,'Left hip joint work','fig');

figure(6);hold on;
errorbar(stiffness,WRhip_pos_ave,WRhip_pos_std);
xlabel('Stiffness(kN/m)');
ylabel('Work(J)');
title('Right hip joint work');
errorbar(stiffness,WRhip_neg_ave,WRhip_neg_std,'r');
hold off
saveas(6,'Right hip joint work','ai');
saveas(6,'Right hip joint work.eps','psc2');
saveas(6,'Right hip joint work','fig');

%% Plot knee work
figure(7);hold on;
errorbar(stiffness,WLkne_pos_ave,WLkne_pos_std);
xlabel('Stiffness(kN/m)');
ylabel('Work(J)');
title('Left knee jointwork');
errorbar(stiffness,WLkne_neg_ave,WLkne_neg_std,'r');
hold off
saveas(7,'Left knee joint work','ai');
saveas(7,'Left knee joint work.eps','psc2');
saveas(7,'Left knee joint work','fig');

figure(8);hold on;
errorbar(stiffness,WRkne_pos_ave,WRkne_pos_std);
xlabel('Stiffness(kN/m)');
ylabel('Work(J)');
title('Right knee work');
errorbar(stiffness,WRkne_neg_ave,WRkne_neg_std,'r');
hold off
saveas(8,'Right knee joint work','ai');
saveas(8,'Right knee joint work.eps','psc2');
saveas(8,'Right knee joint work','fig');

%% Plot ankle work
figure(9);hold on;
errorbar(stiffness,WLank_pos_ave,WLank_pos_std);
xlabel('Stiffness(kN/m)');
ylabel('Work(J)');
title('Left ankle joint work');
errorbar(stiffness,WLank_neg_ave,WLank_neg_std,'r');
hold off
saveas(9,'Left ankle joint work','ai');
saveas(9,'Left ankle joint work','fig');
saveas(9,'Left ankle joint work.eps','psc2');


figure(10);hold on;
errorbar(stiffness,WRank_pos_ave,WRank_pos_std);
xlabel('Stiffness(kN/m)');
ylabel('Work(J)');
title('Right ankle joint work');
errorbar(stiffness,WRank_neg_ave,WRank_neg_std,'r');
hold off
saveas(10,'Right ankle joint work','ai');
saveas(10,'Right ankle joint work.eps','psc2');
saveas(10,'Right ankle joint work','fig');

%% plot hip joint power
figure(11);
plot([0.1:0.1:100]',PRhip_aveR,'LineWidth',2);
xlabel('Time(% per stride)');
ylabel('Power(W)');
legend('6.60 kN/m','18.68 kN/m','38.66 kN/m','58.29 kN/m');
title('Right hip joint power');
saveas(11,'Right hip joint power','ai');
saveas(11,'Right hip joint power.eps','psc2');
saveas(11,'Right hip joint power','fig');

figure(12);
plot([0.1:0.1:100]',PLhip_aveR,'LineWidth',2);
xlabel('Time(% per stride)');
ylabel('Power(W)');
legend('6.60 kN/m','18.68 kN/m','38.66 kN/m','58.29 kN/m');
title('Left hip joint power');
saveas(12,'Left hip joint power','ai');
saveas(12,'Left hip joint power.eps','psc2');
saveas(12,'Left hip joint power','fig');

%% Plot knee joint power
figure(13);
plot([0.1:0.1:100]',PRkne_aveR,'LineWidth',2);
xlabel('Time(% per stride)');
ylabel('Power(W)');
legend('6.60 kN/m','18.68 kN/m','38.66 kN/m','58.29 kN/m');
title('Right knee joint power');
saveas(13,'Right knee joint power','ai');
saveas(13,'Right knee joint power.eps','psc2');
saveas(13,'Right knee joint power','fig');

figure(14);
plot([0.1:0.1:100]',PLkne_aveR,'LineWidth',2);
xlabel('Time(% per stride)');
ylabel('Power(W)');
legend('6.60 kN/m','18.68 kN/m','38.66 kN/m','58.29 kN/m');
title('Left knee joint power');
saveas(14,'Left knee joint power','ai');
saveas(14,'Left knee joint power','fig');
saveas(14,'Left knee joint power.eps','psc2');

%% Plot ankle joint power
figure(15);
plot([0.1:0.1:100]',PRank_aveR,'LineWidth',2);
xlabel('Time(% per stride)');
ylabel('Power(W)');
legend('6.60 kN/m','18.68 kN/m','38.66 kN/m','58.29 kN/m');
title('Right ankle joint power');
saveas(15,'Right ankle joint power','ai');
saveas(15,'Right ankle joint power','fig');
saveas(15,'Right ankle joint power.eps','psc2');

figure(16);
plot([0.1:0.1:100]',PLank_aveR,'LineWidth',2);
xlabel('Time(% per stride)');
ylabel('Power(W)');
legend('6.60 kN/m','18.68 kN/m','38.66 kN/m','58.29 kN/m');
title('Left ankle joint power');
saveas(16,'Left ankle joint power','ai');
saveas(16,'Left ankle joint power','fig');
saveas(16,'Left ankle joint power.eps','psc2');


%% plot hip moment
figure(17);

RRR=[MRhip_aveR(:,1,1) MRhip_aveR(:,1,2) MRhip_aveR(:,1,3) MRhip_aveR(:,1,4) ];
LLL=[MLhip_aveR(:,1,1) MLhip_aveR(:,1,2) MLhip_aveR(:,1,3) MLhip_aveR(:,1,4)];
plot([0.1:0.1:100]',RRR,'LineWidth',2);
xlabel('Time(% per stride)');
ylabel('Moment(Nm)');
legend('6.60 kN/m','18.68 kN/m','38.66 kN/m','58.29 kN/m');
title('Right hip joint moment');
saveas(17,'Right hip joint moment','ai');
saveas(17,'Right hip joint moment.eps','psc2');
saveas(17,'Right hip joint moment','fig');

figure(18);
plot([0.1:0.1:100]',LLL,'LineWidth',2);
xlabel('Time(% per stride)');
ylabel('Moment(Nm)');
legend('6.60 kN/m','18.68 kN/m','38.66 kN/m','58.29 kN/m');
title('Left hip joint moment');
saveas(18,'Left hip joint moment','ai');
saveas(18,'Left hip joint moment','fig');
saveas(18,'Left hip joint moment.eps', 'psc2')
clear RRR LLL
%% plot knee moment
RRR=[MRkne_aveR(:,1,1) MRkne_aveR(:,1,2) MRkne_aveR(:,1,3) MRkne_aveR(:,1,4) ];
LLL=[MLkne_aveR(:,1,1) MLkne_aveR(:,1,2) MLkne_aveR(:,1,3) MLkne_aveR(:,1,4)];

figure(19);
plot([0.1:0.1:100]',RRR,'LineWidth',2);
xlabel('Time(% per stride)');
ylabel('Moment(Nm)');
legend('6.60 kN/m','18.68 kN/m','38.66 kN/m','58.29 kN/m');
title('Right knee joint moment');
saveas(19,'Right knee joint moment','ai');
saveas(19,'Right knee joint moment','fig');
saveas(19,'Right knee joint moment.eps', 'psc2')
figure(20);
plot([0.1:0.1:100]',LLL,'LineWidth',2);
xlabel('Time(% per stride)');
ylabel('Moment(Nm)');
legend('6.60 kN/m','18.68 kN/m','38.66 kN/m','58.29 kN/m');
title('Left knee joint moment');
saveas(20,'Left knee joint moment','ai');
saveas(20,'Left knee joint moment','fig');
saveas(20,'Left knee joint moment.eps', 'psc2')
clear RRR LLL
%% plot ankle moment
RRR=[MRank_aveR(:,1,1) MRank_aveR(:,1,2) MRank_aveR(:,1,3) MRank_aveR(:,1,4) ];
LLL=[MLank_aveR(:,1,1) MLank_aveR(:,1,2) MLank_aveR(:,1,3) MLank_aveR(:,1,4)];

figure(21);
plot([0.1:0.1:100]',RRR,'LineWidth',2);
xlabel('Time(% per stride)');
ylabel('Moment(Nm)');
legend('6.60 kN/m','18.68 kN/m','38.66 kN/m','58.29 kN/m');
title('Right ankle joint moment');
name='Right ankle joint moment';
saveas(21,name,'ai');
saveas(21,name,'fig');
saveas(21,[name '.eps'], 'psc2')
figure(22);
plot([0.1:0.1:100]',LLL,'LineWidth',2);
xlabel('Time(% per stride)');
ylabel('Moment(Nm)');
legend('6.60 kN/m','18.68 kN/m','38.66 kN/m','58.29 kN/m');
title('Left ankle joint moment');
name='Left ankle joint moment';
saveas(22,name,'ai');
saveas(22,name,'fig');
saveas(22,[name '.eps'], 'psc2')

clear RRR LLL
%% plot hip ankle
RRR=[ARhip_aveR(:,1,1) ARhip_aveR(:,1,2) ARhip_aveR(:,1,3) ARhip_aveR(:,1,4) ];
LLL=[ALhip_aveR(:,1,1) ALhip_aveR(:,1,2) ALhip_aveR(:,1,3) ALhip_aveR(:,1,4)];

figure(23);
plot([0.1:0.1:100]',RRR,'LineWidth',2);
xlabel('Time(% per stride)');
ylabel('Degree');
legend('6.60 kN/m','18.68 kN/m','38.66 kN/m','58.29 kN/m');
title('Right hip joint angle');
name='Right hip joint angle';
saveas(23,name,'ai');
saveas(23,name,'fig');
saveas(23,[name '.eps'], 'psc2')

figure(24);
plot([0.1:0.1:100]',LLL,'LineWidth',2);
xlabel('Time(% per stride)');
ylabel('Degree');
legend('6.60 kN/m','18.68 kN/m','38.66 kN/m','58.29 kN/m');
title('Left hip joint angle');
name='Left hip joint angle';
saveas(24,name,'ai');
saveas(24,name,'fig');
saveas(24,[name '.eps'], 'psc2')

clear RRR LLL
%% plot knee angle
RRR=[ARkne_aveR(:,1,1) ARkne_aveR(:,1,2) ARkne_aveR(:,1,3) ARkne_aveR(:,1,4) ];
LLL=[ALkne_aveR(:,1,1) ALkne_aveR(:,1,2) ALkne_aveR(:,1,3) ALkne_aveR(:,1,4)];

figure(25);
plot([0.1:0.1:100]',RRR,'LineWidth',2);
xlabel('Time(% per stride)');
ylabel('Degree');
legend('6.60 kN/m','18.68 kN/m','38.66 kN/m','58.29 kN/m');
title('Right knee joint angle');
name='Right knee joint angle';
saveas(25,name,'ai');
saveas(25,name,'fig');
saveas(25,[name '.eps'], 'psc2')
figure(26);
plot([0.1:0.1:100]',LLL,'LineWidth',2);
xlabel('Time(% per stride)');
ylabel('Degree');
legend('6.60 kN/m','18.68 kN/m','38.66 kN/m','58.29 kN/m');
title('Left knee joint angle');
name='Left knee joint angle';
saveas(26,name,'ai');
saveas(26,name,'fig');
saveas(26,[name '.eps'], 'psc2')
clear RRR LLL
%% plot ankle angle
RRR=[ARank_aveR(:,1,1) ARank_aveR(:,1,2) ARank_aveR(:,1,3) ARank_aveR(:,1,4) ];
LLL=[ALank_aveR(:,1,1) ALank_aveR(:,1,2) ALank_aveR(:,1,3) ALank_aveR(:,1,4)];

figure(27);
plot([0.1:0.1:100]',RRR,'LineWidth',2);
xlabel('Time(% per stride)');
ylabel('Degree');
legend('6.60 kN/m','18.68 kN/m','38.66 kN/m','58.29 kN/m');
title('Right ankle joint angle');
name='Right ankle joint angle';
saveas(27,name,'ai');
saveas(27,[name '.eps'], 'psc2');
saveas(27,name,'fig');

figure(28);
plot([0.1:0.1:100]',LLL,'LineWidth',2);
xlabel('Time(% per stride)');
ylabel('Degree');
legend('6.60 kN/m','18.68 kN/m','38.66 kN/m','58.29 kN/m');
title('Left ankle joint angle');
name='Left ankle joint angle';
saveas(28,name,'ai');
saveas(28,[name '.eps'], 'psc2')
saveas(28,name,'fig');

clear RRR LLL
%% Plot COM work rate
figure(29);
RRR=[RComP_aveR{1} RComP_aveR{2} RComP_aveR{3} RComP_aveR{4}];
LLL=[LComP_aveR{1} LComP_aveR{2} LComP_aveR{3} LComP_aveR{4}];
plot([0.1:0.1:100]',RRR,'LineWidth',2); hold on;
plot([0.1:0.1:100]',LLL,'--','LineWidth',2);
legend('6.60 kN/m','18.68 kN/m','38.66 kN/m','58.29 kN/m');
xlabel('Time(% per stride)');
ylabel('COM work rate(W)');
title('COM work rate');
name='COM work rate';
saveas(29,name,'ai');
saveas(29,[name '.eps'], 'psc2')
saveas(29,name,'fig');

clear RRR LLL

%% Close all
close all



