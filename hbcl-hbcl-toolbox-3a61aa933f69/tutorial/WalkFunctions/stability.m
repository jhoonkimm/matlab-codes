function Z = stability(varargin)
% This function is to search for the step length and width from markers at
% CAL. 
if nargin<2
    error('Too less inputs');
end

RCAL = varargin{1};
LCAL = varargin{2};
if (nargin==2)
trialType = 'Treadmill';
speed = 1.25;
freq = 100;
else 
trialType = varargin{3};
speed = varargin{4};
end

switch (trialType)
    
    case 'Treadmill'
        t=[1/freq:1/freq:size(RCAL,1)/freq];
        RCAL(:,2) = RCAL(:,2)+speed*t';
        LCAL(:,2) = LCAL(:,2)+speed*t';
    case 'Overground'
        t=[1/freq:1/freq:size(RCAL,1)/freq];
end

minRz = min(RCAL(:,3));
minLz = min(LCAL(:,3));




RHS=1;LHS=1;
for i = 2 : size(RCAL,1)-1
    if RCAL(i,3)<RCAL(i-1,3)&&RCAL(i,3)<RCAL(i+1,3)&&RCAL(i,3)<minRz+0.05
        if (RCAL(i,2)-RCAL(RHS(end),2))<0.40
            if RCAL(i,3)<RCAL(RHS(end),3)
                RHS(end)=i;
            end
        else
            RHS=[RHS i];
        end
    end
    if LCAL(i,3)<LCAL(i-1,3)&&LCAL(i,3)<LCAL(i+1,3)&&LCAL(i,3)<minLz+0.05
        if (LCAL(i,2)-LCAL(LHS(end),2))<0.40
            if LCAL(i,3)<LCAL(LHS(end),3)
                LHS(end)=i;
            end
        else
            LHS=[LHS i];
        end
    end
end
RHS=RHS(2:end);
LHS=LHS(2:end);

if RHS(1)>LHS(1)
    LHS=LHS(2:end);
end

if RHS(end)>LHS(end)
    RHS=RHS(1:end-1);
end


figure(111);clf;
subplot(3,1,1);hold on
Rx = line(t,RCAL(:,1),'color','b');
Lx = line(t,LCAL(:,1),'color','r');
line(t(RHS),RCAL(RHS,1),'linestyle','none','marker','o','color','b');
line(t(LHS),LCAL(LHS,1),'linestyle','none','marker','o','color','r');

subplot(3,1,2);hold on
Ry = line(t,RCAL(:,2),'color','b');
Ly = line(t,LCAL(:,2),'color','r');
line(t(RHS),RCAL(RHS,2),'linestyle','none','marker','o','color','b');
line(t(LHS),LCAL(LHS,2),'linestyle','none','marker','o','color','r');

subplot(3,1,3);hold on
Rz = line(t,RCAL(:,3),'color','b');
Lz = line(t,LCAL(:,3),'color','r');
line(t(RHS),RCAL(RHS,3),'linestyle','none','marker','o','color','b');
line(t(LHS),LCAL(LHS,3),'linestyle','none','marker','o','color','r');

for i = 1:length(RHS)-1
    Steps=RCAL(RHS(i):RHS(i+1),:);
    RSTEPS(:,:,i)=interp1([1:size(Steps,1)]',Steps,linspace(1,size(Steps,1),1000));
    RSTEPS(:,:,i)=RSTEPS(:,:,i)-ones(1000,1)*RSTEPS(1,:,i);
end

for i = 1:length(LHS)-1
    Steps=LCAL(LHS(i):LHS(i+1),:);
    LSTEPS(:,:,i)=interp1([1:size(Steps,1)]',Steps,linspace(1,size(Steps,1),1000));
    LSTEPS(:,:,i)=LSTEPS(:,:,i)-ones(1000,1)*LSTEPS(1,:,i);
end

figure(112);clf;
subplot(2,1,1);hold on;
for i = 1:size(RSTEPS,3)
line(RSTEPS(:,2,i),RSTEPS(:,1,i),'color',[0.5 0.5 0.5]);
line(RSTEPS(end,2,i),RSTEPS(end,1,i),'linestyle','none','marker','o','markerfacecolor',[0 0 0],...
    'markeredgecolor','k','markersize',3);
end
axis equal
title('Right')
subplot(2,1,2);hold on;
for i = 1:size(LSTEPS,3)
line(LSTEPS(:,2,i),LSTEPS(:,1,i),'color',[0.5 0.5 0.5]);
line(LSTEPS(end,2,i),LSTEPS(end,1,i),'linestyle','none','marker','o','markerfacecolor',[0 0 0],...
    'markeredgecolor','k','markersize',3);
end
axis equal
title('Left');

for i = 1:length(RHS)-1
    R2LLengths(i,:)=(LCAL(LHS(i),:)-RCAL(RHS(i),:));   
end

for i = 1 :length(LHS)-1
    L2RLengths(i,:)=(-LCAL(LHS(i),:)+RCAL(RHS(i+1),:));
end

figure(113);clf;
subplot(3,1,1);hold on
line([1:size(R2LLengths,1)]',R2LLengths(:,1),'linestyle','none','marker','*','markeredgecolor','b');
line([1:size(L2RLengths,1)]',L2RLengths(:,1),'linestyle','none','marker','*','markeredgecolor','r');
subplot(3,1,2);hold on
line([1:size(R2LLengths,1)]',R2LLengths(:,2),'linestyle','none','marker','*','markeredgecolor','b');
line([1:size(L2RLengths,1)]',L2RLengths(:,2),'linestyle','none','marker','*','markeredgecolor','r');
subplot(3,1,3);hold on
line([1:size(R2LLengths,1)]',R2LLengths(:,3),'linestyle','none','marker','*','markeredgecolor','b');
line([1:size(L2RLengths,1)]',L2RLengths(:,3),'linestyle','none','marker','*','markeredgecolor','r');

Z.RSTEPS=RSTEPS;
Z.LSTEPS=LSTEPS;
Z.R2LLengths=R2LLengths;
Z.L2RLengths=L2RLengths;
% pause
end