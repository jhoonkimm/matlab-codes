function Z = CalC0C(Cop, omega, v, r0)
% this funcion is to calculate the center of curvature of foot for walking. 
% It will find a point with only horizontal velocity 

R=(v(2)+omega(1)*Cop(3))/2/omega(1);

Z.R=R;
r=[Cop(1);Cop(2);R];
Z.r=r;

end