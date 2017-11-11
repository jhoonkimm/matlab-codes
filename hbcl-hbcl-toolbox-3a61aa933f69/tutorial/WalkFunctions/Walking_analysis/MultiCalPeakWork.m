function [W_pos W_neg push_time]=MultiCalPeakWork(P,t,Good,Range)
N=size(P);
if nargin==3
for j=1:N(3)
for i=1:length(Good{j})
    
    [W_pos{j}(i), W_neg{j}(i) ]=CalPeakWork(P(:,Good{j}(i),j),t{j}(i));

end
end
elseif nargin==4
    
for j=1:N(3)
for i=1:length(Good{j})
    [W_pos{j}(i), W_neg{j}(i) push_time{j}(i)]=CalPeakWork_push_time(P(:,Good{j}(i),j),t{j}(i),Range);
i
j
end
end
else
    fprintf('Error: Incorrect number of input arguments');
    return;
end
end