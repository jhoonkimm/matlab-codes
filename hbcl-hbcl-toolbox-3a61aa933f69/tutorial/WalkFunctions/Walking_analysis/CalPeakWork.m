function [W_pos W_neg]=CalPeakWork(Power,t,Range)

if nargin==2
    Range=[1:length(Power)];
elseif nargin==3
    Range=Range;
else
    fprintf('Error: Incorrect number of input arguments');
    return;
end

[C,I]=max(Power(Range));
I=Range(1)+I;
index=I;
sign=1;
while sign>0
    sign=Power(index)*Power(index-1);
    index=index-1;
    Pos_left=index;
end
index=I;
sign=1;
while sign>0
    sign=Power(index)*Power(index+1);
    index=index+1;
    Pos_right=index;
end
index=Pos_left-1;
Neg_right=index;
sign=1;
while sign>0
   sign= Power(index)*Power(index-1);
   Neg_left=index;
   index=index-1;
   if index<=1;
       break
   end
   
   
end
time=t*[1:1000]/1000;
W_pos=trapz(time(Pos_left:Pos_right),Power(Pos_left:Pos_right));
W_neg=trapz(time(Neg_left:Neg_right),Power(Neg_left:Neg_right));


end