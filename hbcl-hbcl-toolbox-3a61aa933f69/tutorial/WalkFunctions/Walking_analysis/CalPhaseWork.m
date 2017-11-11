function [HS RB PL PO]=CalPhaseWork(R,Target,tR)
% This function is to calculate work of target(ex. hip joint power)in 4 phases.
% R is right leg COM work rate started at right heel-strike to next right heel-strike
% Target is the series of power to calculate
% tR is the time period of this strike
% Edit by Paul 2010/10/19


RPO_period=[0 0];
RPL_period=[0 0];
RRB_period=[0 0];
RHS_period=[0 0];
timeR=tR*[1:size(R,1)]/size(R,1);
% % Push-off
% [peak index]=max(R(501:end));
% index=index+500;
% for i=index:1000
%     if R(i)<1e-6
%         Toe_off=i;
%         break
%     end
% end
% for i=index:-1:1
%     if R(i)*R(i+1)<0
%         Push_off=i;
%         break
%     end
% end
% 
% % Rebound
% [peak index]=max(R(1:400));
% for i=index:1000
%     if R(i)*R(i+1)<0
%         Preload=i;
%         break
%     end
% end
% for i=index:-1:1
%     if R(i)*R(i+1)<0
%         Rebound=i;
%         break
%     end
% end
for i=size(R,1):-1:2
   if abs(R(i))>=17&&R(i)>0&&RPO_period(2)==0
       RPO_period(2)=i;
       continue
   end
   if ((R(i)*R(i-1))<0)&&(RPO_period(1)==0)&&(RPO_period(2)~=0)
       RPO_period(1)=i;
       RPL_period(2)=i-1;
       continue
   end
   if ((R(i)*R(i-1))<0)&&RPL_period(1)==0&&RPL_period(2)~=0
       RPL_period(1)=i;
       RRB_period(2)=i-1;
       continue
   end
   if ((R(i)*R(i-1))<0)&&RRB_period(1)==0&&RRB_period(2)~=0
       RRB_period(1)=i;
       RHS_period(2)=i-1;
       RHS_period(1)=1;
       continue
   end
   
end

HS=trapz(timeR(RHS_period(1):RHS_period(2)),Target(RHS_period(1):RHS_period(2)));
RB=trapz(timeR(RRB_period(1):RRB_period(2)),Target(RRB_period(1):RRB_period(2)));
PL=trapz(timeR(RPL_period(1):RPL_period(2)),Target(RPL_period(1):RPL_period(2)));
PO=trapz(timeR(RPO_period(1):RPO_period(2)),Target(RPO_period(1):RPO_period(2)));
% HS=trapz(timeR(1:Rebound),Target(1:Rebound));
% RB=trapz(timeR(Rebound:Preload),Target(Rebound:Preload));
% PL=trapz(timeR(Preload:Push_off),Target(Preload:Push_off));
% PO=trapz(timeR(Push_off:Toe_off),Target(Push_off:Toe_off));
end