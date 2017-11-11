
function [theta cosine_theta]=CalculateVectorAngle(A,B)
% This function calculate the angle between two vectors.
% Edit by Paul 2010/08/26
Length_A=sqrt(A(1)^2+A(2)^2+A(3)^2);
Length_B=sqrt(B(1)^2+B(2)^2+B(3)^2);
cosine_theta=dot(A,B)/Length_A/Length_B;
theta=acos(cosine_theta);


end