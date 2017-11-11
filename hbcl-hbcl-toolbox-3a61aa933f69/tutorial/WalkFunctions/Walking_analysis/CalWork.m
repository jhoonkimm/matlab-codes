function [HS RB PL PO]=CalWork(R,tr)
% Function to calculate for work of four phases from single COM work rate
% series.
% Edit by Paul 2010/10/19
RPO_period=[0 0];
RPL_period=[0 0];
RRB_period=[0 0];
RHS_period=[0 0];

timeR=tr*[1:size(R,1)]/size(R,1);

% [indices] = findPhasesOfGait(R,1)
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

HS=trapz(timeR(RHS_period(1):RHS_period(2)),R(RHS_period(1):RHS_period(2)));
RB=trapz(timeR(RRB_period(1):RRB_period(2)),R(RRB_period(1):RRB_period(2)));
PL=trapz(timeR(RPL_period(1):RPL_period(2)),R(RPL_period(1):RPL_period(2)));
PO=trapz(timeR(RPO_period(1):RPO_period(2)),R(RPO_period(1):RPO_period(2)));
% HS=trapz(timeR(indices(1):indices(2)),R(indices(1):indices(2)));
% RB=trapz(timeR(indices(2):indices(3)),R(indices(2):indices(3)));
% PL=trapz(timeR(indices(3):indices(4)),R(indices(3):indices(4)));
% PO=trapz(timeR(indices(4):indices(5)),R(indices(4):indices(5)));


end
