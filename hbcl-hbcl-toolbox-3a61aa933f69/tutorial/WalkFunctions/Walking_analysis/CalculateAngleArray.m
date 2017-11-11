function[thetas consines]=CalculateAngleArray(A,B)
% This function calculate the angles between two vector arrays
% Edit by Paul 2010/0826
Lengths_A=sqrt(A(:,1).^2+A(:,2).^2+A(:,3).^2);
Lengths_B=sqrt(B(:,1).^2+B(:,2).^2+B(:,3).^2);
consines=dot(A,B,2)./Lengths_A./Lengths_B;
thetas=acos(consines);
end