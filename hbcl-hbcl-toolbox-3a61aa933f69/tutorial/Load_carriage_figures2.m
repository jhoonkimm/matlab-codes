%% Data processing
clear all
close all
clc
% Process individual data
load('result_06232011.mat')
Data=[];
Data_remove_mean=[];
g=9.81;
Ms=[];Lg=[];
normalwalking=[248.46 nan 327.01 313.94 421.94 360.63 437.25 348.05]';
normalLoad = [];
rest_metabolic = [];
%Subject 1
n=1;
token=Amy;
Ms=[Ms token{2}(1)];
Lg=[Lg token{2}(2)];
factor=token{2}(1)*g*sqrt(g*token{2}(2));
Amy_normal=[(token{1}(:,1)+3.8)./token{2}(1) (token{1}(:,2)-Amy_standing)./factor token{1}(:,3)./factor^2];
Data=[Data Amy_normal(:,2)];
normalwalking(n)=normalwalking(n)/factor-Amy_standing/factor;
normalLoad = [normalLoad Amy_normal(:,1)];
rest_metabolic(n) = Amy_standing/factor;
% Subject 2
n=2;
token=Celeste;
Ms=[Ms token{2}(1)];
Lg=[Lg token{2}(2)];
factor=token{2}(1)*g*sqrt(g*token{2}(2));
Celeste_normal=[(token{1}(:,1)+3.8)./token{2}(1) (token{1}(:,2)-Celeste_standing)./factor token{1}(:,3)./factor^2];
Data=[Data Celeste_normal(:,2)];
normalwalking(n)=normalwalking(n)/factor-Celeste_standing/factor;
normalLoad = [normalLoad Celeste_normal(:,1)];
rest_metabolic(n) = Celeste_standing/factor;
% Subject 3
n=3;
token=Huan;
Ms=[Ms token{2}(1)];
Lg=[Lg token{2}(2)];
factor=token{2}(1)*g*sqrt(g*token{2}(2));
Huan_normal=[(token{1}(:,1)+3.8)./token{2}(1) (token{1}(:,2)-Huan_standing)./factor token{1}(:,3)./factor^2];
Data=[Data Huan_normal(:,2)];
normalwalking(n)=normalwalking(n)/factor-Huan_standing/factor;
normalLoad = [normalLoad Huan_normal(:,1)];
rest_metabolic(n) = Huan_standing/factor;
% Subject 4
n=4;
token=Ke;
Ms=[Ms token{2}(1)];
Lg=[Lg token{2}(2)];
factor=token{2}(1)*g*sqrt(g*token{2}(2));
Ke_normal=[(token{1}(:,1)+3.8)./token{2}(1) (token{1}(:,2)-Ke_standing)./factor token{1}(:,3)./factor^2];
Data=[Data Ke_normal(:,2)];
normalwalking(n)=normalwalking(n)/factor-Ke_standing/factor;
normalLoad = [normalLoad Ke_normal(:,1)];
rest_metabolic(n) = Ke_standing/factor;
% Subject 5
n=5;
token=Paul;
Ms=[Ms token{2}(1)];
Lg=[Lg token{2}(2)];
factor=token{2}(1)*g*sqrt(g*token{2}(2));
Paul_normal=[(token{1}(:,1)+3.8)./token{2}(1) (token{1}(:,2)-Paul_standing)./factor token{1}(:,3)./factor^2];
Data=[Data Paul_normal(:,2)];

normalwalking(n)=normalwalking(n)/factor-Paul_standing/factor;

normalLoad = [normalLoad Paul_normal(:,1)];
rest_metabolic(n) = Paul_standing/factor;
% Subject 6
n=6;
token=Shioa;
Ms=[Ms token{2}(1)];
Lg=[Lg token{2}(2)];
factor=token{2}(1)*g*sqrt(g*token{2}(2));
Shioa_normal=[(token{1}(:,1)+3.8)./token{2}(1) (token{1}(:,2)-Shioa_standing)./factor token{1}(:,3)./factor^2];
Data=[Data Shioa_normal(:,2)];
normalwalking(n)=normalwalking(n)/factor-Shioa_standing/factor;
normalLoad = [normalLoad Shioa_normal(:,1)];
rest_metabolic(n) = Shioa_standing/factor;
% Subject 7
n=7;
token=Watson;
Ms=[Ms token{2}(1)];
Lg=[Lg token{2}(2)];
factor=token{2}(1)*g*sqrt(g*token{2}(2));
Watson_normal=[(token{1}(:,1)+3.8)./token{2}(1) (token{1}(:,2)-Watson_standing)./factor token{1}(:,3)./factor^2];
Data=[Data Watson_normal(:,2)];
normalwalking(n)=normalwalking(n)/factor-Watson_standing/factor;
normalLoad = [normalLoad Watson_normal(:,1)];
rest_metabolic(n) = Watson_standing/factor;
% Subject 8
n=8;
token=Mike;
Ms=[Ms token{2}(1)];
Lg=[Lg token{2}(2)];
factor=token{2}(1)*g*sqrt(g*token{2}(2));
Mike_normal=[(token{1}(:,1)+3.8)./token{2}(1) (token{1}(:,2)-Mike_standing)./factor token{1}(:,3)./factor^2];
Data=[Data [Mike_normal([1 2],2);nan;Mike_normal([3],2)]];
normalwalking(n)=normalwalking(n)/factor-Watson_standing/factor;
normalLoad = [normalLoad [Mike_normal([1 2],1);nan;Mike_normal([3],1)]];
rest_metabolic(n) = Mike_standing/factor;
load('LoadCarriageExperimentResult_include_huan_v4.mat')
%%
figure(1);clf;
width = 600; height =200;
set(1,'position',[500 500 width height]);
data_numbers ={};
data_numbers{1,1} = ['8'];
data_numbers{2,1} = [num2str(mean(Ms),'%5.3f') '+\-' num2str(std(Ms),'%5.3f')]; 
data_numbers{3,1} = [num2str(mean(Lg),'%5.3f') '+\-' num2str(std(Lg),'%5.3f')]; 
data_numbers{4,1} = [num2str(mean(rest_metabolic),'%5.3f') '+\-' num2str(std(rest_metabolic),'%5.3f')];
data_numbers{5,1} = [num2str(mean(normalwalking([1 3:8])),'%5.3f') '+\-' num2str(std(normalwalking([1 3:8])),'%5.3f')];
data_numbers{6,1} = [num2str(mean(steplength_normal),'%5.3f') '+\-' num2str(std(steplength_normal),'%5.3f')];
step_frequency = DS_time_normal./DS_percentage_normal*100;
data_numbers{7,1} = [num2str(mean(step_frequency),'%5.3f') '+\-' num2str(std(step_frequency),'%5.3f')];

table1 = uitable('position',[0 0 width 200],'Data',data_numbers, 'ColumnName', {'mean +/- s.d.'},...
                        'Rowname',{'Subject #','Body mass (kg)','Leg length (m)','Rest metabolic rate','No load walking metabolic rate','No load step length','No load step frequency (1/s)'});


%% Fig 2: Ground reaction force and COM work rate

Totalmass = [ones(size(normalLoad,1),1) normalLoad+1];
figure(2);clf;
set(2,'position',[100 100 1000 400],...
    'name','Ground reaction force and COM work rate');
subplot(1,2,2);hold on;
ax1=gca;
set(ax1,'ycolor',[1 1 1],'ytick',[],'xtick',[0:20:100],'xlim',[0 100],...
    'box','off','color','none');
xlabel('Time(% stride)','fontsize',14);
axis_position = get(ax1,'position');
ax2=axes('Position',axis_position+[0 axis_position(4)*0.05 0 0] ,...
         'YAxisLocation','left',...
         'XAxisLocation','bot',...
         'xlim',[0 100],'xtick',[],'xcolor',[1 1 1],...
         'Color','none',...
         'ylim',[-0.15 0.15]);
line([0.1:0.1:100]',RComP_normal_aveR,'linewidth',2,'color',[190 190 190 ]./255);
line([0.1:0.1:100]',RComP_aveR,'linewidth',2); 
line([0.1:0.1:100]',LComP_aveR,'linewidth',2,'linestyle','--');
line([0.1:0.1:100]',LComP_normal_aveR,'linewidth',2,'color',[190 190 190 ]./255,'linestyle','--');
line([0 100],[0 0],'color','k');
l1 = legend('Normal','15 lb','25 lb','35 lb','45 lb' );
set(get(ax2,'ylabel'),'string','W','fontsize',14)
set(l1,'location','southeast')
ylabel('COM work rate','fontsize',14)
range=[-0.15 0.15]*Mass_ave*g*sqrt(leglength_ave*g);
ax3=axes('Position',axis_position+[0 axis_position(4)*0.05 0 0] ,...
         'YAxisLocation','right',...
         'XAxisLocation','bot',...
         'xlim',[0 100],'xtick',[],'xcolor',[1 1 1],...
         'Color','none',...
         'ylim',range);
ylabel('W','fontsize',14);

subplot(1,2,1);hold on;

line([0.1:0.1:100]',RGRFz_normal_ave,'linewidth',2,'color',[190 190 190 ]./255);
line([0.1:0.1:100]',RGRFz_ave,'linewidth',2); hold on;
line([0.1:0.1:100]',LGRFz_normal_ave,'linewidth',2,'color',[190 190 190 ]./255);
line([0.1:0.1:100]',LGRFz_ave,'linewidth',2,'linestyle','--');
xlabel('Time(% stride)','fontsize',14);
ylabel('Vertical GRF','fontsize',14);
hold off
ax1=gca;
axis_position = get(ax1,'position');
set(ax1,'color','none','position',axis_position+[0 axis_position(4)*0.05 0 0])
range=get(ax1,'ylim')*Mass_ave*g;
l1=legend('Normal','15 lb','25 lb','35 lb','45 lb' );

ax2=axes('Position',axis_position+[0 axis_position(4)*0.05 0 0],...
         'YAxisLocation','right',...
         'XAxisLocation','top',...
         'xtick',[],...
         'Color','none',...
         'ylim',range);
set(get(ax2,'ylabel'),'string','N','fontsize',14)
set(l1,'location','southeast')

%% Fig 3: 3 x 3
figure(3);clf
set(3,'position',[100 100 1500 750],...
    'name','3 by 3')
subplot(3,3,7);hold on
ax1=gca;
set(ax1,'ycolor',[1 1 1],'ytick',[],'xtick',[0:20:100],'xlim',[0 100],...
    'box','off','color','none');
xlabel('Time(% stride)','fontsize',14);
axis_position = get(ax1,'position');
ax2=axes('Position',axis_position+[0 axis_position(4)*0.1 0 0] ,...
         'YAxisLocation','left',...
         'XAxisLocation','bot',...
         'xlim',[0 100],'xtick',[],'xcolor',[1 1 1],...
         'Color','none',...
         'ylim',[-0.1 0.25]);
line([0.1:0.1:100]',PRhip_aveR,'linewidth',linewidth);
legend('15 lb','25 lb','35 lb','45 lb' );
line([0.1:0.1:100]',PRhip_normal_aveR,'linewidth',linewidth,'color',[190 190 190 ]./255);
ylabel('Power','fontsize',14);
range = get(ax2,'ylim');
range=range*Mass_ave*g*sqrt(leglength_ave*g);
ax3=axes('Position',get(ax2,'position'),...
         'YAxisLocation','right',...
         'xtick',[],'xcolor',[1 1 1],...
         'Color','none',...
         'ylim',range);
set(get(ax3,'ylabel'),'string','W','fontsize',10)

clear ax1 ax2
hold off

subplot(3,3,8);hold on
ax1=gca;
set(ax1,'ycolor',[1 1 1],'ytick',[],'xtick',[0:20:100],'xlim',[0 100],...
    'box','off','color','none');
xlabel('Time(% stride)','fontsize',14);
axis_position = get(ax1,'position');
ax2=axes('Position',axis_position+[0 axis_position(4)*0.1 0 0] ,...
         'YAxisLocation','left',...
         'XAxisLocation','bot',...
         'xlim',[0 100],'xtick',[],'xcolor',[1 1 1],...
         'Color','none',...
         'ylim',[-0.1 0.25]);
line([0.1:0.1:100]',PRkne_aveR,'linewidth',linewidth);
axis([0 100 -0.1 0.25]);
line([0.1:0.1:100]',PRkne_normal_aveR,'linewidth',linewidth,'color',[190 190 190 ]./255);
ylabel('Power','fontsize',14);
range = get(ax2,'ylim');
range=range*Mass_ave*g*sqrt(leglength_ave*g);
ax3=axes('Position',get(ax2,'position'),...
         'YAxisLocation','right',...
         'xtick',[],'xcolor',[1 1 1],...
         'Color','none',...
         'ylim',range);
set(get(ax3,'ylabel'),'string','W','fontsize',10)

clear ax1 ax2
hold off

subplot(3,3,9);hold on
ax1=gca;
set(ax1,'ycolor',[1 1 1],'ytick',[],'xtick',[0:20:100],'xlim',[0 100],...
    'box','off','color','none');
xlabel('Time(% stride)','fontsize',14);
axis_position = get(ax1,'position');
ax2=axes('Position',axis_position+[0 axis_position(4)*0.1 0 0] ,...
         'YAxisLocation','left',...
         'XAxisLocation','bot',...
         'xlim',[0 100],'xtick',[],'xcolor',[1 1 1],...
         'Color','none',...
         'ylim',[-0.1 0.25]);
line([0.1:0.1:100]',PRank_aveR,'linewidth',linewidth);
axis([0 100 -0.1 0.25]);
line([0.1:0.1:100]',PRank_normal_aveR,'linewidth',linewidth,'color',[190 190 190 ]./255);
ylabel('Power','fontsize',14);
range = get(ax2,'ylim');
range=range*Mass_ave*g*sqrt(leglength_ave*g);
ax3=axes('Position',get(ax2,'position'),...
         'YAxisLocation','right',...
         'xtick',[],'xcolor',[1 1 1],...
         'Color','none',...
         'ylim',range);
set(get(ax3,'ylabel'),'string','W','fontsize',10)

clear ax1 ax2
hold off



% plot hip moment
figure(3);subplot(3,3,4);hold on
line([0.1:0.1:100]',-MRhip_normal_aveR,'linewidth',linewidth,'color',[190 190 190 ]./255);
line([0.1:0.1:100]',-MRhip_aveR,'linewidth',linewidth);
axis([0 100 -0.1 0.25]);
ylabel('Moment','fontsize',14);
ax1=gca;
range=get(ax1,'ylim')*Mass_ave*g*leglength_ave;
set(ax1,'color','none','xcolor',[1 1 1]);
ax2=axes('Position',get(ax1,'position'),...
         'YAxisLocation','right',...
         'xtick',[],'xcolor',[1 1 1],...
         'Color','none',...
         'ylim',range);
set(get(ax2,'ylabel'),'string','N-m','fontsize',10)
hold off

figure(3);subplot(3,3,5);hold on
line([0.1:0.1:100]',MRkne_aveR,'linewidth',linewidth);
 line([0.1:0.1:100]',-MRkne_normal_aveR,'linewidth',linewidth,'color',[190 190 190 ]./255);
axis([0 100 -0.1 0.25])
ylabel('','fontsize',14);
ax1=gca;
set(ax1,'color','none','xcolor',[1 1 1]);
range=get(ax1,'ylim')*Mass_ave*g*leglength_ave;
ax2=axes('Position',get(ax1,'position'),...
         'YAxisLocation','right',...
         'xtick',[],'xcolor',[1 1 1],...
         'Color','none',...
         'ylim',range);
set(get(ax2,'ylabel'),'string','N-m','fontsize',10)

hold off

figure(3);subplot(3,3,6);hold on
 line([0.1:0.1:100]',-MRank_normal_aveR,'linewidth',linewidth,'color',[190 190 190 ]./255);

line([0.1:0.1:100]',-MRank_aveR,'linewidth',linewidth);
axis([0 100 -0.1 0.25])
name='Right ankle joint moment';
ax1=gca;
set(ax1,'color','none','xcolor',[1 1 1]);
range=get(ax1,'ylim')*Mass_ave*g*leglength_ave;
ax2=axes('Position',get(ax1,'position'),...
         'YAxisLocation','right',...
         'xtick',[],'xcolor',[1 1 1],...
         'Color','none',...
         'ylim',range);
set(get(ax2,'ylabel'),'string','N-m','fontsize',10)
hold off

%plot hip ankle

subplot(3,3,1);hold on
 line([0.1:0.1:100]',-(ARhip_normal_aveR-ARhip_normal_aveR(1,1)),'linewidth',linewidth,'color',[190 190 190 ]./255);
line([0.1:0.1:100]',-(ARhip_aveR-ARhip_aveR(1,1)),'linewidth',linewidth);
axis([0 100 -90 60]);
ylabel('Angle (degree)','fontsize',14);
title('Hip','fontsize',14);
name='Right hip joint angle';
ax1=gca;
set(ax1,'color','none','xcolor',[1 1 1],'ytick',[-90:30:60]);


subplot(3,3,2);hold on
 line([0.1:0.1:100]',-(ARkne_normal_aveR-ARkne_normal_aveR(1,1)),'linewidth',linewidth,'color',[190 190 190 ]./255);
line([0.1:0.1:100]',(ARkne_aveR-ARkne_aveR(1,1)),'linewidth',linewidth);
axis([0 100 -90 60]);
title('Knee','fontsize',14);
ax1=gca;
range=get(ax1,'ylim')*180/pi;
set(ax1,'color','none','xcolor',[1 1 1],'ytick',[-90:30:60]);

% plot ankle angle
hold off
subplot(3,3,3);hold on
line([0.1:0.1:100]',-(ARank_normal_aveR-ARank_normal_aveR(1,1)),'linewidth',linewidth,'color',[190 190 190 ]./255);
line([0.1:0.1:100]',-(ARank_aveR-ARank_aveR(1,1)),'linewidth',linewidth);
axis([0 100 -90 60]);
title('Ankle','fontsize',14);
ax1=gca;
set(ax1,'color','none','xcolor',[1 1 1],'ytick',[-90:30:60]);
hold off





%% Fig 4: Work per step
figure(4);clf;hold all
set(4,'position',[100 100 1000 375],...
    'name','Work per step');
subplot(1,2,1);hold all;
axis([0.9 1.5 -0.1 0.1])
ax1=gca;
set(ax1,'ycolor',[1 1 1],'ytick',[],'xtick',[0.9:0.1:1.5],...
    'box','off','color','none');
xlabel('Total mass','fontsize',14);
axis_position = get(ax1,'position');
ax2=axes('Position',axis_position+[0 axis_position(4)*0.05 0 0] ,...
         'YAxisLocation','left',...
         'XAxisLocation','bot',...
         'xlim',[0.9 1.5],'xtick',[],'xcolor',[1 1 1],...
         'Color','none',...
         'ylim',[-0.15 0.15]);
hold on
ylabel('Work/step','fontsize',14)

x_data = [normalLoad+1];
y_data = PO_ave_normal+RB_ave_normal;

T=4;
[b_load_vs_com_pos,bint_load_vs_com_pos,r_load_vs_com_pos,rint_load_vs_com_pos,stats_load_vs_com_pos] =...
    regress(reshape(y_data',[],1),[reshape(x_data',[],1)  blkdiag(ones(T,1), ones(T,1) ,ones(T,1), ones(T,1) ,ones(T,1) ,ones(T,1), ones(T,1),ones(T,1))]);
COM_positive=plot(x_data',(y_data-b_load_vs_com_pos(2:end)*ones(1,4))'+mean(b_load_vs_com_pos(2:end)),'cd','markerfacecolor','c');
plot(ones(1,8),PO_normalwalking_normal+RB_normalwalking_normal-b_load_vs_com_pos(2:end)'+mean(b_load_vs_com_pos(2:end)),'cd')
line([1 1.5],[1 1.5]*b_load_vs_com_pos(1)+mean(b_load_vs_com_pos(2:end)),'linewidth',2,'color','c')

y_data = Bio_work_pos;

T=4;
[b_load_vs_bio_pos,bint_load_vs_bio_pos,r_load_vs_com_pos,rint_load_vs_com_pos,stats_load_vs_bio_pos] =...
    regress(reshape(y_data',[],1),[reshape(x_data',[],1)  blkdiag(ones(T,1), ones(T,1) ,ones(T,1), ones(T,1) ,ones(T,1) ,ones(T,1), ones(T,1),ones(T,1))]);
Bio=plot(x_data',(y_data-b_load_vs_bio_pos(2:end)*ones(1,4))'+mean(b_load_vs_bio_pos(2:end)),'rs','markerfacecolor','r');
plot(ones(1,8),Bio_work_pos_normal_dimensionless -b_load_vs_bio_pos(2:end)'+mean(b_load_vs_bio_pos(2:end)),'rs');
line([1 1.5],[1 1.5]*b_load_vs_bio_pos(1)+mean(b_load_vs_bio_pos(2:end)),'linewidth',2,'color','r');

y_data = Bio_work_neg;

T=4;
[b_load_vs_bio_neg,bint_load_vs_bio_neg,r_load_vs_com_pos,rint_load_vs_com_pos,stats_load_vs_bio_neg] =...
    regress(reshape(y_data',[],1),[reshape(x_data',[],1)  blkdiag(ones(T,1), ones(T,1) ,ones(T,1), ones(T,1) ,ones(T,1) ,ones(T,1), ones(T,1),ones(T,1))]);
Bio=plot(x_data',(y_data-b_load_vs_bio_neg(2:end)*ones(1,4))'+mean(b_load_vs_bio_neg(2:end)),'rs','markerfacecolor','r');
plot(ones(1,8),Bio_work_neg_normal_dimensionless -b_load_vs_bio_neg(2:end)'+mean(b_load_vs_bio_neg(2:end)),'rs');
line([1 1.5],[1 1.5]*b_load_vs_bio_neg(1)+mean(b_load_vs_bio_neg(2:end)),'linewidth',2,'color','r')

% MiddleAxis(0.0025,[0.9:0.1:1.4],[0.9 1.5],0.015,-0.0075)
plot([0.9 1.5],[0 0],'k')
rangey=get(ax2,'ylim')*Mass_ave*g;
ax3=axes('Position',get(ax2,'position'),...
         'YAxisLocation','right',...
         'Color','none',...
         'ylim',rangey, 'ytick',[-100 :50:100],...
         'xlim',[0.9 1.5],'xtick',[],'xcolor',[1 1 1]);
ylabel('J','fontsize',14)
l1=legend([COM_positive(1) Bio(1)],'COM work','Positive of sum of joint work');
set(l1,'location','best')
%
subplot(1,2,2);hold all
axis([0.9 1.5 -0.1 0.1])
ax1=gca;
set(ax1,'ycolor',[1 1 1],'ytick',[],'xtick',[0.9:0.1:1.5],...
    'box','off','color','none');
xlabel('Total mass','fontsize',14);
axis_position = get(ax1,'position');
ax2=axes('Position',axis_position+[0 axis_position(4)*0.05 0 0] ,...
         'YAxisLocation','left',...
         'XAxisLocation','bot',...
         'xlim',[0.9 1.5],'xtick',[],'xcolor',[1 1 1],...
         'box','off',...
         'Color','none',...
         'ylim',[-0.1 0.1],'ytick',[-0.1:0.05:0.1]);

hold all;
ylabel('Work/step','fontsize',14)

% y_data = Ank_work_pos+Kne_work_pos+Hip_work_pos;
% T = 4;
% [b_load_vs_sum_of_pos,bint_load_vs_sum_of_pos,r_load_vs_sum_of_pos,rint_load_vs_sum_of_pos,stats_load_vs_sum_of_pos] =...
%     regress(reshape(y_data',[],1),[reshape(x_data',[],1)  blkdiag(ones(T,1), ones(T,1) ,ones(T,1), ones(T,1) ,ones(T,1) ,ones(T,1), ones(T,1),ones(T,1))]);
% 
% Sum_of_joints=plot(x_data',(y_data-b_load_vs_sum_of_pos(2:end)*ones(1,4))'+mean(b_load_vs_sum_of_pos(2:end)),'ko','markerfacecolor','k');
% 
% plot(ones(1,8), Ank_work_pos_normal+Kne_work_pos_normal+Hip_work_pos_normal,'ko');
% line([1 1.5],[1 1.5]*b_load_vs_sum_of_pos(1)+mean(b_load_vs_sum_of_pos(2:end)),'linewidth',2,'color','k');
% 
% y_data = Ank_work_neg+Hip_work_neg+Kne_work_neg;
% T = 4;
% [b_load_vs_sum_of_neg,bint_load_vs_sum_of_neg,r_load_vs_sum_of_neg,rint_load_vs_sum_of_neg,stats_load_vs_sum_of_neg] =...
%     regress(reshape(y_data',[],1),[reshape(x_data',[],1)  blkdiag(ones(T,1), ones(T,1) ,ones(T,1), ones(T,1) ,ones(T,1) ,ones(T,1), ones(T,1),ones(T,1))]);
% 
% Sum_of_joints=plot(x_data',(y_data-b_load_vs_sum_of_neg(2:end)*ones(1,4))'+mean(b_load_vs_sum_of_neg(2:end)),'ko','markerfacecolor','k');
% 
% plot(ones(1,8), Ank_work_neg_normal+Hip_work_neg_normal+Kne_work_neg_normal,'ko');
% line([1 1.5],[1 1.5]*b_load_vs_sum_of_neg(1)+mean(b_load_vs_sum_of_neg(2:end)),'linewidth',2,'color','k')
% 

y_data = Hip_work_pos;
T = 4;
[b_hip_work_vs_load,bint_hip_work_vs_load,r_hip_work_vs_load,rint_hip_work_vs_load,stats_hip_work_vs_load] =...
    regress(reshape(y_data',[],1),[reshape(x_data',[],1)  blkdiag(ones(T,1), ones(T,1) ,ones(T,1), ones(T,1) ,ones(T,1) ,ones(T,1), ones(T,1),ones(T,1))]);

plot(x_data',(y_data-b_hip_work_vs_load(2:end)*ones(1,4))'+mean(b_hip_work_vs_load(2:end)),'^','color',[0 0.5 0],'markerfacecolor',[0 0.5 0])

plot(ones(1,8), Hip_work_pos_normal_dimensionless -b_hip_work_vs_load(2:end)'+mean(b_hip_work_vs_load(2:end)),'^');
line([1 1.5],[1 1.5]*b_hip_work_vs_load(1)+mean(b_hip_work_vs_load(2:end)),'linewidth',2,'color',[0 0.5 0])

y_data = Hip_work_neg;
T = 4;
[b_hip_work_vs_load_neg,bint_hip_work_vs_load_neg,r_hip_work_vs_load_neg,rint_hip_work_vs_load_neg,stats_hip_work_vs_load_neg] =...
    regress(reshape(y_data',[],1),[reshape(x_data',[],1)  blkdiag(ones(T,1), ones(T,1) ,ones(T,1), ones(T,1) ,ones(T,1) ,ones(T,1), ones(T,1),ones(T,1))]);

Hip = plot(x_data',(y_data-b_hip_work_vs_load_neg(2:end)*ones(1,4))'+mean(b_hip_work_vs_load_neg(2:end)),'^','color',[0 0.5 0],'markerfacecolor',[0 0.5 0]);

plot(ones(1,8), Hip_work_neg_normal_dimensionless -b_hip_work_vs_load_neg(2:end)'+mean(b_hip_work_vs_load_neg(2:end)),'^','color',[0 0.5 0]);
line([1 1.5],[1 1.5]*b_hip_work_vs_load_neg(1)+mean(b_hip_work_vs_load_neg(2:end)),'linewidth',2,'color',[0 0.5 0]);

y_data = Kne_work_pos;
T = 4;
[b_kne_work_vs_load,bint_kne_work_vs_load,r_kne_work_vs_load,rint_kne_work_vs_load,stats_kne_work_vs_load] =...
    regress(reshape(y_data',[],1),[reshape(x_data',[],1)  blkdiag(ones(T,1), ones(T,1) ,ones(T,1), ones(T,1) ,ones(T,1) ,ones(T,1), ones(T,1),ones(T,1))]);

plot(x_data',(y_data-b_kne_work_vs_load(2:end)*ones(1,4))'+mean(b_kne_work_vs_load(2:end)),'d','color','m','markerfacecolor','r');

plot(ones(1,8), Kne_work_pos_normal_dimensionless -b_kne_work_vs_load(2:end)'+mean(b_kne_work_vs_load(2:end)),'md');
line([1 1.5],[1 1.5]*b_kne_work_vs_load(1)+mean(b_kne_work_vs_load(2:end)),'linewidth',2,'color','m');

y_data = Kne_work_neg;
T = 4;
[b_kne_work_vs_load_neg,bint_kne_work_vs_load_neg,r_kne_work_vs_load,rint_kne_work_vs_load,stats_kne_work_vs_load_neg] =...
    regress(reshape(y_data',[],1),[reshape(x_data',[],1)  blkdiag(ones(T,1), ones(T,1) ,ones(T,1), ones(T,1) ,ones(T,1) ,ones(T,1), ones(T,1),ones(T,1))]);

Knee = plot(x_data',(y_data-b_kne_work_vs_load_neg(2:end)*ones(1,4))'+mean(b_kne_work_vs_load_neg(2:end)),'d','color','m','markerfacecolor','r');

plot(ones(1,8), Kne_work_neg_normal-b_kne_work_vs_load_neg(2:end)'+mean(b_kne_work_vs_load_neg(2:end)),'md');
line([1 1.5],[1 1.5]*b_kne_work_vs_load_neg(1)+mean(b_kne_work_vs_load_neg(2:end)),'linewidth',2,'color','m');

y_data = Ank_work_pos;
T = 4;
[b_ank_work_vs_load,bint_ank_work_vs_load,r_ank_workrate_vs_load,rint_ank_workrate_vs_load,stats_ank_work_vs_load] =...
    regress(reshape(y_data',[],1),[reshape(x_data',[],1)  blkdiag(ones(T,1), ones(T,1) ,ones(T,1), ones(T,1) ,ones(T,1) ,ones(T,1), ones(T,1),ones(T,1))]);

plot(x_data',(y_data-b_ank_work_vs_load(2:end)*ones(1,4))'+mean(b_ank_work_vs_load(2:end)),'s','color','b','markerfacecolor','b')

plot(ones(1,8), Ank_work_pos_normal_dimensionless -b_ank_work_vs_load(2:end)'+mean(b_ank_work_vs_load(2:end)),'bs');
line([1 1.5],[1 1.5]*b_ank_work_vs_load(1)+mean(b_ank_work_vs_load(2:end)),'linewidth',2,'color','b')

y_data = Ank_work_neg;
T = 4;
[b_ank_work_vs_load_neg,bint_ank_work_vs_load_neg,r_ank_workrate_vs_load,rint_ank_workrate_vs_load,stats_ank_work_vs_load_neg] =...
    regress(reshape(y_data',[],1),[reshape(x_data',[],1)  blkdiag(ones(T,1), ones(T,1) ,ones(T,1), ones(T,1) ,ones(T,1) ,ones(T,1), ones(T,1),ones(T,1))]);

Ankle = plot(x_data',(y_data-b_ank_work_vs_load_neg(2:end)*ones(1,4))'+mean(b_ank_work_vs_load_neg(2:end)),'s','color','b','markerfacecolor','b');

plot(ones(1,8), Ank_work_neg_normal_dimensionless -b_ank_work_vs_load_neg(2:end)'+mean(b_ank_work_vs_load_neg(2:end)),'bs');
line([1 1.5],[1 1.5]*b_ank_work_vs_load_neg(1)+mean(b_ank_work_vs_load_neg(2:end)),'linewidth',2,'color','b');

plot([0.9 1.5],[0 0],'k')

rangey=get(ax2,'ylim')*Mass_ave*g;
ax3=axes('Position',get(ax2,'position'),...
         'YAxisLocation','right',...
         'Color','none',...
         'ylim',rangey, 'ytick',[-100 :50:100],...
         'xlim',[0.9 1.5],'xtick',[],'xcolor',[1 1 1]);
ylabel('J','fontsize',14)
l1 = legend([ Hip(1) Knee(1) Ankle(1)],'Hip work','Knee work','Ankle work');
set(l1,'location','Best')

%
figure(41);clf
width = 600; height =200;
set(41,'position',[500 500 width height]);
data_numbers ={};
data_numbers{1,1} = [ num2str(b_load_vs_com_pos(1),'%4.3f') '+/-'  num2str(bint_load_vs_com_pos(1,2)-b_load_vs_com_pos(1),'%4.3f') ];
data_numbers{1,2} = [ num2str(b_load_vs_bio_pos(1),'%4.3f') '+/-'  num2str(bint_load_vs_bio_pos(1,2)-b_load_vs_bio_pos(1),'%4.3f') ];
data_numbers{1,3} = [ num2str(b_load_vs_bio_neg(1),'%4.3f') '+/-'  num2str(bint_load_vs_bio_neg(1,2)-b_load_vs_bio_neg(1),'%4.3f') ];
data_numbers{2,1} = [ num2str(mean(b_load_vs_com_pos(2:end)),'%4.3f') '+/-'  num2str(std(b_load_vs_com_pos(2:end)),'%4.3f') ];
data_numbers{2,2} = [ num2str(mean(b_load_vs_bio_pos(2:end)),'%4.3f') '+/-'  num2str(std(b_load_vs_bio_pos(2:end)),'%4.3f') ];
data_numbers{2,3} = [ num2str(mean(b_load_vs_bio_neg(2:end)),'%4.3f') '+/-'  num2str(std(b_load_vs_bio_neg(2:end)),'%4.3f') ];
data_numbers{3,1} = num2str(stats_load_vs_com_pos(1),'%4.3f');
data_numbers{3,2} = num2str(stats_load_vs_bio_pos(1),'%4.3f');
data_numbers{3,3} = num2str(stats_load_vs_bio_neg(1),'%4.3f');
data_numbers{4,1} = num2str(stats_load_vs_com_pos(3),'%4.3e');
data_numbers{4,2} = num2str(stats_load_vs_bio_pos(3),'%4.3e');
data_numbers{4,3} = num2str(stats_load_vs_bio_neg(3),'%4.3e');
table1 = uitable('position',[0 height-100 width 100],'Data',data_numbers, 'ColumnName', {'*COM work', '*Pos of sum', '*Neg of sum'},...
                        'Rowname',{'Slope','Offset','R^2','p value'});
%                 
data_numbers ={};
data_numbers{1,1} = [ num2str(b_hip_work_vs_load(1),'%4.3f') '+/-'  num2str(bint_hip_work_vs_load(1,2)-b_hip_work_vs_load(1),'%4.3f') ];
data_numbers{1,2} = [ num2str(b_hip_work_vs_load_neg(1),'%4.3f') '+/-'  num2str(bint_hip_work_vs_load_neg(1,2)-b_hip_work_vs_load_neg(1),'%4.3f') ];
data_numbers{1,3} = [ num2str(b_kne_work_vs_load(1),'%4.3f') '+/-'  num2str(bint_kne_work_vs_load(1,2)-b_kne_work_vs_load(1),'%4.3f') ];
data_numbers{1,4} = [ num2str(b_kne_work_vs_load_neg(1),'%4.3f') '+/-'  num2str(bint_kne_work_vs_load_neg(1,2)-b_kne_work_vs_load_neg(1),'%4.3f') ];
data_numbers{1,5} = [ num2str(b_ank_work_vs_load(1),'%4.3f') '+/-'  num2str(bint_ank_work_vs_load(1,2)-b_ank_work_vs_load(1),'%4.3f') ];
data_numbers{1,6} = [ num2str(b_ank_work_vs_load_neg(1),'%4.3f') '+/-'  num2str(bint_ank_work_vs_load_neg(1,2)-b_ank_work_vs_load_neg(1),'%4.3f') ];


data_numbers{2,1} = [ num2str(mean(b_hip_work_vs_load(2:end)),'%4.3f') '+/-'  num2str(std(b_hip_work_vs_load(2:end)),'%4.3f') ];
data_numbers{2,2} = [ num2str(mean(b_hip_work_vs_load_neg(2:end)),'%4.3f') '+/-'  num2str(std(b_hip_work_vs_load_neg(2:end)),'%4.3f') ];
data_numbers{2,3} = [ num2str(mean(b_kne_work_vs_load(2:end)),'%4.3f') '+/-'  num2str(std(b_kne_work_vs_load(2:end)),'%4.3f') ];
data_numbers{2,4} = [ num2str(mean(b_kne_work_vs_load_neg(2:end)),'%4.3f') '+/-'  num2str(std(b_kne_work_vs_load_neg(2:end)),'%4.3f') ];
data_numbers{2,5} = [ num2str(mean(b_ank_work_vs_load(2:end)),'%4.3f') '+/-'  num2str(std(b_ank_work_vs_load(2:end)),'%4.3f') ];
data_numbers{2,6} = [ num2str(mean(b_ank_work_vs_load_neg(2:end)),'%4.3f') '+/-'  num2str(std(b_ank_work_vs_load_neg(2:end)),'%4.3f') ];



data_numbers{3,1} = num2str(stats_hip_work_vs_load(1),'%4.3f');
data_numbers{3,2} = num2str(stats_hip_work_vs_load_neg(1),'%4.3f');
data_numbers{3,3} = num2str(stats_kne_work_vs_load(1),'%4.3f');
data_numbers{3,4} = num2str(stats_kne_work_vs_load_neg(1),'%4.3f');
data_numbers{3,5} = num2str(stats_ank_work_vs_load(1),'%4.3f');
data_numbers{3,6} = num2str(stats_ank_work_vs_load_neg(1),'%4.3f');

data_numbers{4,1} = num2str(stats_hip_work_vs_load(3),'%4.3e');
data_numbers{4,2} = num2str(stats_hip_work_vs_load_neg(3),'%4.3e');
data_numbers{4,3} = num2str(stats_kne_work_vs_load(3),'%4.3e');
data_numbers{4,4} = num2str(stats_kne_work_vs_load_neg(3),'%4.3e');
data_numbers{4,5} = num2str(stats_ank_work_vs_load(3),'%4.3e');
data_numbers{4,6} = num2str(stats_ank_work_vs_load_neg(3),'%4.3e');


table2 = uitable('position',[0 height-200 width 100],'Data',data_numbers, 'ColumnName', {'Hip +','*Hip -','*Knee +','Knee -','*Ankle +','*Ankle -'},...
                        'Rowname',{'Slope','Offset','R^2','p value'});
                                  
figure(42);clf
width = 700; height =200;
set(42,'position',[500 500 width height],'name','Normal walking data: work/step');

data_numbers ={};
data_numbers{1,1} =[num2str(mean(PO_normalwalking_normal+RB_normalwalking_normal),'%4.3f') '+/-' num2str(std(PO_normalwalking_normal+RB_normalwalking_normal),'%4.3f')];
data_numbers{1,2} =[num2str(mean(Bio_work_pos_normal_dimensionless ),'%4.3f') '+/-' num2str(std(Bio_work_pos_normal_dimensionless ),'%4.3f')];
data_numbers{1,3} =[num2str(mean(Bio_work_neg_normal_dimensionless ),'%4.3f') '+/-' num2str(std(Bio_work_neg_normal_dimensionless ),'%4.3f')];

data_numbers{2,1} =[num2str(mean(PO_normalwalking+RB_normalwalking),'%4.3f') '+/-' num2str(std(PO_normalwalking+RB_normalwalking),'%4.3f')];
data_numbers{2,2} =[num2str(mean(Bio_work_pos_normal),'%4.3f') '+/-' num2str(std(Bio_work_pos_normal),'%4.3f')];
data_numbers{2,3} =[num2str(mean(Bio_work_neg_normal),'%4.3f') '+/-' num2str(std(Bio_work_neg_normal),'%4.3f')];


table1 = uitable('position',[0 height-100 width 100],'Data',data_numbers, 'ColumnName', {'*COM work', '*Pos of sum', '*Neg of sum'},...
                   'Rowname',{'dimensionless','J'},'ColumnWidth',{90});
data_numbers ={};
Ank_work_pos_normal_dimensionless
data_numbers{1,1} =[num2str(mean(Hip_work_pos_normal_dimensionless),'%4.3f') '+/-' num2str(std(Hip_work_pos_normal_dimensionless),'%4.3f')];
data_numbers{1,2} =[num2str(mean(Hip_work_neg_normal_dimensionless ),'%4.3f') '+/-' num2str(std(Hip_work_neg_normal_dimensionless ),'%4.3f')];
data_numbers{1,3} =[num2str(mean(Kne_work_pos_normal_dimensionless ),'%4.3f') '+/-' num2str(std(Kne_work_pos_normal_dimensionless ),'%4.3f')];
data_numbers{1,4} =[num2str(mean(Kne_work_neg_normal_dimensionless ),'%4.3f') '+/-' num2str(std(Kne_work_neg_normal_dimensionless ),'%4.3f')];
data_numbers{1,5} =[num2str(mean(Ank_work_pos_normal_dimensionless ),'%4.3f') '+/-' num2str(std(Ank_work_pos_normal_dimensionless ),'%4.3f')];
data_numbers{1,6} =[num2str(mean(Ank_work_neg_normal_dimensionless ),'%4.3f') '+/-' num2str(std(Ank_work_neg_normal_dimensionless ),'%4.3f')];

data_numbers{2,1} =[num2str(mean(Hip_work_pos_normal),'%4.3f') '+/-' num2str(std(Hip_work_pos_normal),'%4.3f')];
data_numbers{2,2} =[num2str(mean(Hip_work_neg_normal ),'%4.3f') '+/-' num2str(std(Hip_work_neg_normal ),'%4.3f')];
data_numbers{2,3} =[num2str(mean(Kne_work_pos_normal ),'%4.3f') '+/-' num2str(std(Kne_work_pos_normal ),'%4.3f')];
data_numbers{2,4} =[num2str(mean(Kne_work_neg_normal ),'%4.3f') '+/-' num2str(std(Kne_work_neg_normal ),'%4.3f')];
data_numbers{2,5} =[num2str(mean(Ank_work_pos_normal ),'%4.3f') '+/-' num2str(std(Ank_work_pos_normal ),'%4.3f')];
data_numbers{2,6} =[num2str(mean(Ank_work_neg_normal ),'%4.3f') '+/-' num2str(std(Ank_work_neg_normal ),'%4.3f')];

table3 = uitable('position',[0 height-200 width 100],'Data',data_numbers, 'ColumnName', {'Hip +','*Hip -','*Knee +','Knee -','*Ankle +','*Ankle -'},...
                       'Rowname',{'dimensionless','J'},'ColumnWidth',{90});

%% Fig 5: Rate of mechanical work (Comparison of normalization)
figure(5);clf;hold all;

set(5,'position',[100 100 1000 375],...
    'name','Rate of mechanical work (Comparison of normalization)')
subplot(1,2,1);hold all
T=4;
y_data = PO_rate+RB_rate;
[b_load_vs_comrate_pos,bint_load_vs_comrate_pos,r_load_vs_com_pos,rint_load_vs_com_pos,stats_load_vs_comrate_pos] =...
    regress(reshape(y_data',[],1),[reshape(x_data',[],1)  blkdiag(ones(T,1), ones(T,1) ,ones(T,1), ones(T,1) ,ones(T,1) ,ones(T,1), ones(T,1),ones(T,1))]);
COM_positive=plot(x_data',(y_data-b_load_vs_comrate_pos(2:end)*ones(1,4))'+mean(b_load_vs_comrate_pos(2:end)),'.','markersize',20);
plot(ones(1,8),PO_workrate_normal+RB_workrate_normal-b_load_vs_comrate_pos(2:end)'+mean(b_load_vs_comrate_pos(2:end)),'ko');
line([1 1.5],[1 1.5]*b_load_vs_comrate_pos(1)+mean(b_load_vs_comrate_pos(2:end)),'linewidth',2,'color','k')
axis([0.9 1.5 0 0.05])
ylabel('Rate of mechanical work','fontsize',14)
ax1=gca;
set(ax1,'ytick',[0:0.01:0.05],'xtick',[])
rangey=get(ax1,'ylim')*Mass_ave*g*sqrt(leglength_ave*g);
ax2=axes('Position',get(ax1,'position'),...
         'YAxisLocation','right',...
         'Color','none',...
         'ylim',rangey,...
         'xlim',[0.9 1.5],...
         'xtick',[0.9:0.1:1.5],...
         'ytick',[0:25:100]);
set(get(ax2,'ylabel'),'string','W','fontsize',14)
xlabel('Total mass','fontsize',14);
subplot(1,2,2);hold all;
y_data = PO_rate_W+RB_rate_W;
[b_load_vs_com_pos_W,bint_load_vs_com_pos_W,r_load_vs_com_pos,rint_load_vs_com_pos,stats_load_vs_com_pos_W] =...
    regress(reshape(y_data',[],1),[reshape(x_data',[],1) ones(size(reshape(x_data',[],1)))]);
COM_positive=plot(x_data',y_data','.','markersize',20);
plot(ones(1,8),PO_workrate_normal_W+RB_workrate_normal_W,'ko');
line([1 1.5],[1 1.5]*b_load_vs_com_pos_W(1)+mean(b_load_vs_com_pos_W(2:end)),'linewidth',2,'color','k')
axis([0.9 1.5 0 0.05*Mass_ave*g*sqrt(leglength_ave*g)])
ylabel('Rate of mechanical work (W)','fontsize',14)
ax1=gca;
set(ax1,'ytick',[0:25:100],'xtick',[0.9:0.1:1.5])

ax2=axes('Position',get(ax1,'position'),...
         'YAxisLocation','right',...
         'Color','none',...
         'ylim',[0 1],...
         'xlim',[0.9 1.5],...
         'xtick',[0.9:0.1:1.5],...
         'ytick',[]);
xlabel('Total mass','fontsize',14);


figure(51);clf
width = 600; height =100;
set(51,'position',[500 500 width height]);
data_numbers ={};
data_numbers{1,1} = [ num2str(b_load_vs_comrate_pos(1),'%4.3f') '+/-'  num2str(bint_load_vs_comrate_pos(1,2)-b_load_vs_comrate_pos(1),'%4.3f') ];
data_numbers{1,2} = [ num2str(b_load_vs_com_pos_W(1),'%4.3f') '+/-'  num2str(bint_load_vs_com_pos_W(1,2)-b_load_vs_com_pos_W(1),'%4.3f') ];

data_numbers{2,1} = [ num2str(mean(b_load_vs_comrate_pos(2:end)),'%4.3f') '+/-'  num2str(std(b_load_vs_comrate_pos(2:end)),'%4.3f') ];
data_numbers{2,2} = [ num2str(mean(b_load_vs_com_pos_W(2:end)),'%4.3f') '+/-'  num2str(std(b_load_vs_com_pos_W(2:end)),'%4.3f') ];

data_numbers{3,1} = num2str(stats_load_vs_comrate_pos(1),'%4.3f');
data_numbers{3,2} = num2str(stats_load_vs_com_pos_W(1),'%4.3f');

data_numbers{4,1} = num2str(stats_load_vs_comrate_pos(3),'%4.3e');
data_numbers{4,2} = num2str(stats_load_vs_com_pos_W(3),'%4.3e');

table1 = uitable('position',[0 height-100 width 100],'Data',data_numbers, 'ColumnName', {' *Non-dimensionalized', ' Regular'     },...
                        'Rowname',{'Slope','Offset','R^2','p value'});
                    
figure(52);clf
width = 600; height =200;
set(52,'position',[500 500 width height],'name','Normal walking data: work rate');

data_numbers ={};
data_numbers{1,1} =[num2str(mean(PO_workrate_normal+RB_workrate_normal),'%4.3f') '+/-' num2str(std(PO_workrate_normal+RB_workrate_normal),'%4.3f')];

data_numbers{2,1} =[num2str(mean(PO_workrate_normal_W+PO_workrate_normal_W),'%4.3f') '+/-' num2str(std(PO_workrate_normal_W+PO_workrate_normal_W),'%4.3f')];


table1 = uitable('position',[0 height-100 width 100],'Data',data_numbers, 'ColumnName', {'COM work rate'},...
                   'Rowname',{'dimensionless','J'},'ColumnWidth',{90});

%% Fig 6: Work per step vs Mv^2 s^2
load('MaxData.mat')
figure(6);clf;hold all
set(6,'name','Work per step vs Mv^2 s^2')
y_data = squeeze(condavgwrk2(:,17,:)+condavgwrk2(:,19,:))./mglrep;
x_data = (condavgspds./sqrtglrep).^2 .* (condavgstplngs./leglenrep).^2;
Old = plot(x_data,y_data,'.','markersize',20,'color',[0.8 0.8 0.8]);
v_normal = 1.25./sqrt(g*ls)'*ones(1,4);
x_data = (normalLoad+1).*v_normal.^2.*steplength.^2;
y_data = PO_ave_normal+RB_ave_normal;
[b_com_pos_vs_mvs,bint_com_pos_vs_mvs,r,rint,stats_com_pos_vs_mvs] =...
    regress(reshape(y_data',[],1),[reshape(x_data',[],1)  blkdiag(ones(T,1), ones(T,1) ,ones(T,1), ones(T,1) ,ones(T,1) ,ones(T,1), ones(T,1),ones(T,1))]);
New = plot(x_data',(y_data-b_com_pos_vs_mvs(2:end)*ones(1,4))'+mean(b_com_pos_vs_mvs(2:end)),'k.','markersize',20);
plot(v_normal(:,1)'.^2.*steplength_normal.^2, PO_normalwalking_normal+RB_normalwalking_normal,'ko');
line([0 0.18],[0 0.18]*b_com_pos_vs_mvs(1)+mean(b_com_pos_vs_mvs(2:end)),'linewidth',2,'color','k');
l1 = legend([New(1) Old(1)],'Vary mass','Vary speed and step length');
set(l1,'location','best') ;
ylabel('Work/step','fontsize',14)
xlabel('M V^2 s^2','fontsize',14)

figure(61);clf;hold all
set(61,'name','Work per step vs Mv^2 s^2')
v_normal = 1.25./sqrt(g*ls)'*ones(1,4);
x_data = (normalLoad+1).*v_normal.^2.*steplength.^2;
y_data = PO_ave_normal+RB_ave_normal;
[b_com_pos_vs_mvs,bint,r,rint,stats_com_pos_vs_mvs] =...
    regress(reshape(y_data',[],1),[reshape(x_data',[],1)  blkdiag(ones(T,1), ones(T,1) ,ones(T,1), ones(T,1) ,ones(T,1) ,ones(T,1), ones(T,1),ones(T,1))]);
New = plot(x_data',(y_data-b_com_pos_vs_mvs(2:end)*ones(1,4))'+mean(b_com_pos_vs_mvs(2:end)),'.','markersize',20);
plot(v_normal(:,1)'.^2.*steplength_normal.^2, PO_normalwalking_normal+RB_normalwalking_normal-b_com_pos_vs_mvs(2:end)'+mean(b_com_pos_vs_mvs(2:end)),'ko');
line([0.05 0.14],[0.05 0.14]*b_com_pos_vs_mvs(1)+mean(b_com_pos_vs_mvs(2:end)),'linewidth',2,'color','k');
% l1 = legend([New(1) Old(1)],'Vary mass','Vary speed and step length');
set(l1,'location','best') ;
ylabel('Work/step','fontsize',14)
xlabel('M V^2 s^2','fontsize',14)
axis([0 0.15 0 0.15])
ax1 = gca;
axis_position = get(ax1,'position');

set(ax1,'xlim',[0 0.15], 'xtick',[0:0.05:0.15],...
        'ylim',[0 0.15],'ytick',[0:0.05:0.15],...
         'position',axis_position-[0 0 axis_position(3)*0.1 axis_position(4)*0.1],...
         'Color','none');
rangex = get(ax1,'xlim')*Mass_ave*leglength_ave^3*g;
rangey = get(ax1,'ylim')*Mass_ave*g*leglength_ave;
ax2=axes('Position',get(ax1,'position'),...
         'YAxisLocation','right',...
         'Xaxislocation','top',...
         'Color','none',...
         'ylim',rangey,...
         'xlim',rangex,...
         'xtick',[0:20:100],...
         'ytick',[0:25:100]);
xlabel('kg m^4 s^{-2}','fontsize',14)
ylabel('J','fontsize',14)

figure(62);clf
width = 600; height =100;
set(62,'position',[500 500 width height]);
data_numbers ={};
data_numbers{1,1} = [ num2str(b_com_pos_vs_mvs(1),'%4.3f') '+/-'  num2str(bint_com_pos_vs_mvs(1,2)-b_com_pos_vs_mvs(1),'%4.3f') ];

data_numbers{2,1} = [ num2str(mean(b_com_pos_vs_mvs(2:end)),'%4.3f') '+/-'  num2str(std(b_com_pos_vs_mvs(2:end)),'%4.3f') ];

data_numbers{3,1} = num2str(stats_com_pos_vs_mvs(1),'%4.3f');

data_numbers{4,1} = num2str(stats_com_pos_vs_mvs(3),'%4.3e');

table1 = uitable('position',[0 height-100 width 100],'Data',data_numbers, 'ColumnName', {' *COM work /step vs MV^2s^2'     },...
                        'Rowname',{'Slope','Offset','R^2','p value'});

%% Fig 7 : Net mechanical rate vs rate of mechanical work
figure(7);clf;hold all 
set(7,'position',[100 100 1000 375],...
    'name','Net mechanical rate vs rate of mechanical work')
subplot(1,2,2);hold all
% y_data = [normalwalking(1:8) Data'];
% x_data = [ones(8,1) normalLoad+1];
% [b_meta_vs_load,bint_meta_vs_load,r_meta_vs_load,rint_meta_vs_load,stats_meta_vs_load] =...
%     regress(reshape(y_data',[],1),[reshape(x_data',[],1)  blkdiag(ones(5,1), ones(5,1) ,ones(5,1), ones(5,1) ,ones(5,1) ,ones(5,1), ones(5,1),ones(5,1))]);
y_data = [Data'];
x_data = [normalLoad+1];
[b_meta_vs_load,bint_meta_vs_load,r_meta_vs_load,rint_meta_vs_load,stats_meta_vs_load] =...
    regress(reshape(y_data',[],1),[reshape(x_data',[],1)  blkdiag(ones(4,1), ones(4,1) ,ones(4,1), ones(4,1) ,ones(4,1) ,ones(4,1),ones(4,1), ones(4,1))]);
line([1 1.5],[1 1.5]*b_meta_vs_load(1)+mean(b_meta_vs_load(2:end)),'linewidth',2,'color','k');
plot(x_data',(y_data-b_meta_vs_load(2:end)*ones(1,4))'+mean(b_meta_vs_load(2:end)),'.','markersize',20)

plot(ones(8,1),normalwalking(1:8)-b_meta_vs_load(2:end)+ mean(b_meta_vs_load(2:end)),'ko')

M_ave=mean(Ms);
Lg_ave=mean(Lg);
ylabel('Net metabolic rate','fontsize',14);
xlabel('Total mass','fontsize',14);

  axis([0.9 1.5 0 0.25])
ax1=gca;
axis_position = get(ax1,'position');

set(ax1,'xtick',[0.9:0.1:1.5],...
    'Color','none',...
    'position',axis_position-[0 0 0 axis_position(4)*0.1]);
range=get(ax1,'ylim')*M_ave*g*sqrt(g*Lg_ave);
ax2=axes('Position',get(ax1,'position'),...
         'YAxisLocation','right',...
         'XAxisLocation','top',...
         'xlim',[0.9 1.5],'xtick',[],...
         'Color','none',...
         'ylim',range,'ytick',[0:100:500]);
set(get(ax2,'ylabel'),'string','W','fontsize',14)

hold off
subplot(1,2,1);hold all
x_data = Data';
y_data = PO_rate+RB_rate;
T=4;
[b_meta_vs_com_pos,bint_meta_vs_com_pos,r_meta_vs_com_pos,rint_meta_vs_com_pos,stats_meta_vs_com_pos] =...
    regress(reshape(y_data',[],1),[reshape(x_data',[],1)  blkdiag(ones(T,1), ones(T,1) ,ones(T,1), ones(T,1) ,ones(T,1) ,ones(T,1), ones(T,1),ones(T,1))]);
COM_positive=plot(x_data',(y_data-b_meta_vs_com_pos(2:end)*ones(1,4))'+mean(b_meta_vs_com_pos(2:end)),'.','markersize',20);
plot(normalwalking,PO_workrate_normal+RB_workrate_normal-b_meta_vs_com_pos(2:end)'+mean(b_meta_vs_com_pos(2:end)),'ko')
line([0.01,0.25],[0.01 0.25]*b_meta_vs_com_pos(1)+mean(b_meta_vs_com_pos(2:end)),'linewidth',2,'color','k')
axis([0 0.25 0 0.05])
ylabel('Rate of mechanical work','fontsize',14);
xlabel('Net metabolic rate','fontsize',14);

ax1=gca;
axis_position = get(ax1,'position');
set(ax1,'ytick',[0:0.01:0.05],...
    'Color','none',...
    'position',axis_position-[0 0 0 axis_position(4)*0.1]);
range=get(ax1,'ylim')*M_ave*g*sqrt(g*Lg_ave);
rangex = get(ax1,'xlim')*M_ave*g*sqrt(g*Lg_ave); 
ax2=axes('Position',get(ax1,'position'),...
         'YAxisLocation','right',...
         'XAxisLocation','top',...
         'xlim',rangex,'ytick',[0:25:100],...
         'Color','none',...
         'ylim',range,'xtick',[0:100:500]);
set(get(ax2,'ylabel'),'string','W','fontsize',14)
set(get(ax2,'xlabel'),'string','W','fontsize',14)

figure(71);clf
width = 600; height =100;
set(71,'position',[500 500 width height]);
data_numbers ={};
data_numbers{1,1} = [ num2str(b_meta_vs_com_pos(1),'%4.3f') '+/-'  num2str(bint_meta_vs_com_pos(1,2)-b_meta_vs_com_pos(1),'%4.3f') ];
data_numbers{1,2} = [ num2str(b_meta_vs_load(1),'%4.3f') '+/-'  num2str(bint_meta_vs_load(1,2)-b_meta_vs_load(1),'%4.3f') ];

data_numbers{2,1} = [ num2str(mean(b_meta_vs_com_pos(2:end)),'%4.3f') '+/-'  num2str(std(b_meta_vs_com_pos(2:end)),'%4.3f') ];
data_numbers{2,2} = [ num2str(mean(b_meta_vs_load(2:end)),'%4.3f') '+/-'  num2str(std(b_meta_vs_load(2:end)),'%4.3f') ];

data_numbers{3,1} = num2str(stats_meta_vs_com_pos(1),'%4.3f');
data_numbers{3,2} = num2str(stats_meta_vs_load(1),'%4.3f');

data_numbers{4,1} = num2str(stats_meta_vs_com_pos(3),'%4.3e');
data_numbers{4,2} = num2str(stats_meta_vs_load(3),'%4.3e');

table1 = uitable('position',[0 height-100 width 100],'Data',data_numbers, 'ColumnName', {'*COM work rate vs metabolic','*Metabolic vs load'     },...
                        'Rowname',{'Slope','Offset','R^2','p value'});


%% Fig 8: Step length, step time and double support time
figure(8);clf;hold all
set(8,'position',[100 100 1000 375],...
    'name','Step length, step time and double support time');
subplot(1,2,1);hold all
x_data = [normalLoad+1];
y_data =steplength;

T=4;
[b_steplength_vs_load,bint_steplength_vs_load,r,rint,stats_steplength_vs_load] =...
    regress(reshape(y_data',[],1),[reshape(x_data',[],1)  blkdiag(ones(T,1), ones(T,1) ,ones(T,1), ones(T,1) ,ones(T,1) ,ones(T,1), ones(T,1),ones(T,1))]);
plot(x_data',(y_data-b_steplength_vs_load(2:end)*ones(1,4))'+mean(b_steplength_vs_load(2:end)),'.','markersize',20)
plot(ones(1,8),steplength_normal-b_steplength_vs_load(2:end)'+mean(b_steplength_vs_load(2:end)),'ko')
%  plot([ones(8,1) x_data]',[steplength_normal' y_data]','.-','markersize',20)
%  plot(ones(1,8),steplength_normal,'ko')

line([1 1.4],[1 1.4]*b_steplength_vs_load(1)+mean(b_steplength_vs_load(2:end)),'linewidth',2,'color','k')

axis([0.9 1.5 0.5 0.8])  
ylabel('Step length','fontsize',14) 
xlabel('Total mass','fontsize',14)
ax1 = gca;
set(ax1,'Color','none')
range = get(ax1,'ylim')*leglength_ave;
ax2 = axes('position',get(ax1,'position'),...
    'yaxislocation','right',...
    'ylim',range,'ytick',[0.5:0.1:0.8],...
    'xcolor',[1 1 1],...
    'color','none',...
    'xtick',[]);
ylabel('M','fontsize',14)



subplot(1,2,2);hold all
y_data =DS_time;

T=4;
[b_DS_time_vs_load,bint_DS_time_vs_load,r,rint,stats_DS_time_vs_load] =...
    regress(reshape(y_data',[],1),[reshape(x_data',[],1)  blkdiag(ones(T,1), ones(T,1) ,ones(T,1), ones(T,1) ,ones(T,1) ,ones(T,1), ones(T,1),ones(T,1))]);
DS_index = plot(x_data',(y_data-b_DS_time_vs_load(2:end)*ones(1,4))'+mean(b_DS_time_vs_load(2:end)),'b^','markerfacecolor','b');
plot(ones(1,8),DS_time_normal-b_DS_time_vs_load(2:end)'+mean(b_DS_time_vs_load(2:end)),'b^')
line([1 1.4],[1 1.4]*b_DS_time_vs_load(1)+mean(b_DS_time_vs_load(2:end)),'linewidth',2,'color','b')

y_data =DS_time./DS_percentage*100;

T=4;
[b_step_time_vs_load,bint_step_time_vs_load,r,rint,stats_step_time_vs_load] =...
    regress(reshape(y_data',[],1),[reshape(x_data',[],1)  blkdiag(ones(T,1), ones(T,1) ,ones(T,1), ones(T,1) ,ones(T,1) ,ones(T,1), ones(T,1),ones(T,1))]);
ST_index = plot(x_data',(y_data-b_step_time_vs_load(2:end)*ones(1,4))'+mean(b_step_time_vs_load(2:end)),'kd','markerfacecolor','k');
plot(ones(1,8),DS_time_normal./DS_percentage_normal*100-b_step_time_vs_load(2:end)'+mean(b_step_time_vs_load(2:end)),'kd');
line([1 1.4],[1 1.4]*b_step_time_vs_load(1)+mean(b_step_time_vs_load(2:end)),'linewidth',2,'color','k');

axis([0.9 1.5 0 0.6])  
ylabel('Time(s)','fontsize',14) 
xlabel('Total mass','fontsize',14)
set(gca,'color','none')
l1 = legend([ST_index(1) DS_index(1)],'Step time','Double support time');
set(l1,'location','best');

figure(81);clf
width = 600; height =100;
set(81,'position',[500 500 width height]);
data_numbers ={};
data_numbers{1,1} = [ num2str(b_steplength_vs_load(1),'%4.3f') '+/-'  num2str(bint_steplength_vs_load(1,2)-b_steplength_vs_load(1),'%4.3f') ];
data_numbers{1,2} = [ num2str(b_step_time_vs_load(1),'%4.3f') '+/-'  num2str(bint_step_time_vs_load(1,2)-b_step_time_vs_load(1),'%4.3f') ];
data_numbers{1,3} = [ num2str(b_DS_time_vs_load(1),'%4.3f') '+/-'  num2str(bint_DS_time_vs_load(1,2)-b_DS_time_vs_load(1),'%4.3f') ];

data_numbers{2,1} = [ num2str(mean(b_steplength_vs_load(2:end)),'%4.3f') '+/-'  num2str(std(b_steplength_vs_load(2:end)),'%4.3f') ];
data_numbers{2,2} = [ num2str(mean(b_step_time_vs_load(2:end)),'%4.3f') '+/-'  num2str(std(b_step_time_vs_load(2:end)),'%4.3f') ];
data_numbers{2,3} = [ num2str(mean(b_DS_time_vs_load(2:end)),'%4.3f') '+/-'  num2str(std(b_DS_time_vs_load(2:end)),'%4.3f') ];

data_numbers{3,1} = num2str(stats_steplength_vs_load(1),'%4.3f');
data_numbers{3,2} = num2str(stats_step_time_vs_load(1),'%4.3f');
data_numbers{3,3} = num2str(stats_DS_time_vs_load(1),'%4.3f');

data_numbers{4,1} = num2str(stats_steplength_vs_load(3),'%4.3e');
data_numbers{4,2} = num2str(stats_step_time_vs_load(3),'%4.3e');
data_numbers{4,3} = num2str(stats_DS_time_vs_load(3),'%4.3e');

table1 = uitable('position',[0 height-100 width 100],'Data',data_numbers, 'ColumnName', {'Step length','Step time','*Double support time'     },...
                        'Rowname',{'Slope','Offset','R^2','p value'});
figure(82);clf
width = 600; height =200;
set(82,'position',[500 500 width height],'name','Normal walking data: work rate');

data_numbers ={};
data_numbers{1,1} =[num2str(mean(steplength_normal),'%4.3f') '+/-' num2str(std(steplength_normal),'%4.3f')];

data_numbers{1,2} =[num2str(mean(DS_time_normal./DS_percentage_normal*100),'%4.3f') '+/-' num2str(std(DS_time_normal./DS_percentage_normal*100),'%4.3f')];
data_numbers{1,3} =[num2str(mean(DS_time_normal),'%4.3f') '+/-' num2str(std(DS_time_normal),'%4.3f')];
data_numbers{2,1} =[num2str(mean(steplength_normal_meter),'%4.3f') '+/-' num2str(std(steplength_normal_meter),'%4.3f')];


table1 = uitable('position',[0 height-100 width 100],'Data',data_numbers, 'ColumnName', {'Step length','Step time','DS time'},...
                   'Rowname',{'dimensionless','M'},'ColumnWidth',{90});

%%
%%

