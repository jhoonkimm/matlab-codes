

%[RComP LComP RComP_ave LComP_ave RComP_std LComP_std velocity]=CalComPower (RGRFxr, RGRFyr, RGRFzr,LGRFxl,LGRFyl,LGRFzl,t,Vel);
M=(RGRFzr+LGRFzl)/9.81;
M=mean(mean(mean(M,1),2),3);
clear vx vy vz
vx=zeros(size(velocity,1),size(velocity,3));
vy=vx;vz=vx;
g=9.81;
for i =1:size(velocity,3)
    vx(:,i)=velocity(:,1,i);
    vy(:,i)=velocity(:,2,i);
    vz(:,i)=velocity(:,3,i);
end
RComP_ave_norm=RComP_ave./M/g/Vel;
LComP_ave_norm=LComP_ave./M/g/Vel;
plot(RComP_ave_norm,'LineWidth',2)
hold on
plot(LComP_ave_norm,'--','LineWidth',2)
legend('1/8','3/16','1/4','3/8')
title('Varun')
figure;
plot(vy,'LineWidth',2);
legend('1/8','3/16','1/4','3/8')
title('Varun y-velocity')
figure;plot(vz,'LineWidth',2);
legend('1/8','3/16','1/4','3/8')
title('Varun z-velocity')

figure;plot(lgrfz_ave,'--','LineWidth',2);
hold on;
plot(rgrfz_ave,'LineWidth',2);
hold off
legend('1/8','3/16','1/4','3/8')
title('Varun z-grf')

figure;plot(lgrfy_ave,'--','LineWidth',2);
hold on;
plot(rgrfy_ave,'LineWidth',2);
hold off
legend('1/8','3/16','1/4','3/8')
title('Varun y-grf')


