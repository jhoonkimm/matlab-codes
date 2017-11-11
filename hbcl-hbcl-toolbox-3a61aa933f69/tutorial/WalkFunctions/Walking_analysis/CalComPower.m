function [RComP LComP RComP_ave LComP_ave RComP_std LComP_std velocity]=CalComPower (Rx, Ry, Rz,Lx,Ly,Lz,t,Good,Vel,varargin)
StepNumber=size(Good,2);
fileNum=size(Rx,3);
%% Optional input
opt_val=varargin;
while length(opt_val)>=2
    opt=opt_val{1};
    val=opt_val{2};
    
    switch opt
        case 'StepNumber'
            StepNumber=val;
       
    end
    if length(opt_val)==2
        break
    else
        opt_val=opt_val{3:end};
    end
end
RComP=zeros(size(Rx,1),StepNumber,size(Rx,3));
LComP=zeros(size(Lx,1),StepNumber,size(Lx,3));

RComP_ave=zeros(size(Rx,1),size(Rx,3));
LComP_ave=zeros(size(Lx,1),size(Lx,3));

RComP_std=zeros(size(Rx,1),size(Rx,3));
LComP_std=zeros(size(Lx,1),size(Lx,3));

for k=1:fileNum
  for j=1:size(Good,2)
    RGRF=[Rx(:,Good(j),k) Ry(:,Good(j),k) Rz(:,Good(j),k)];
    LGRF=[Lx(:,Good(j),k) Ly(:,Good(j),k) Lz(:,Good(j),k)];
    [RComP(:,j,k) ,LComP(:,j,k), v_temp(:,:,j)]=calculateComWorkRateForSteadyStateGait(RGRF,LGRF,Vel,t(Good(j),k),'subtractVelocitySlope',1);
  end
  
  RComP_ave(:,k)=mean(RComP,2);
  RComP_std(:,k)=std(RComP,0,2);
  velocity(:,:,k)=mean(v_temp,3);
  clear v_temp;
  
  LComP_ave(:,k)=mean(LComP,2);
  LComP_std(:,k)=std(LComP,0,2);
end


end