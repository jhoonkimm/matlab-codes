function [HS_ave HS_std HSs RB_ave RB_std RBs PL_ave PL_std PLs PO_ave PO_std POs]=MultiCalWork(RComP,Target,tr,GoodR,N)
% MultiCalWork calculates Collision, Rebound, Preload and Push-off work
% from COM work rate. RComP and LComP represent COM work for two legs which
% are (numbers of data per step) by (numbers of steps) by (numbers of
% trials). RComP starts at right heel-strike to right heel-strike. LComP 
% starts at left heel-strike to left heel-strike. tr and tl represent 
% stride time which are (numbers of steps) by (numbers of trial).
% Update by Paul 2010/10/19

for i =1:N
    for j=1:25%length(GoodR{i})
        i
        j
        [HSs{i}(j) RBs{i}(j) PLs{i}(j) POs{i}(j)]=CalPhaseWork(RComP{i}(:,GoodR{i}(j)),Target{i}(:,GoodR{i}(j)),tr{i}(GoodR{i}(j)));
    end

%     for j=1:length(GoodL{i})
%         
%         [HSs{1}{i}(j) RBs{1}{i}(j) PLs{1}{i}(j) POs{1}{i}(j)]=CalWork(LComP{i}(:,GoodL{i}(j)),tl{i}(GoodL{i}(j)));
%     end
    HS_ave(i)=mean(HSs{i} );
    PL_ave(i)=mean(PLs{i} );
    RB_ave(i)=mean(RBs{i});
    PO_ave(i)=mean(POs{i} );
    HS_std(i)=std(HSs{i} );
    PL_std(i)=std(PLs{i} );
    RB_std(i)=std(RBs{i} );
    PO_std(i)=std(POs{i} );

end
end