function d = diffshc(x)

% DIFFSHC same as diff but output vector is same length as input
% diffshc(x) returns a vector d representing the
% difference between sequential values of x, similar
% to diff.m, but of the same length as x.  The missing
% value is approximated using polynomial approximation.
% This function is useful for discrete derivative 
% approximation.
% 
% In order to prevent a time shift, the differece is
% first shifted into the middle through a linear 
% interpolation, then both ends are addended using
% a 2nd order polynomial approximation.  Therefore,
% the minimum number of rows required is 5.
% 
% Created by Steve Collins, shc@umich.edu
% Updated by Karl Zelik 11/19/09


if size(x,2) > size(x,1)
    x=x';
end

d0 = diff(x);
d = zeros(size(x));

for n = 1:size(x,2)
    d(2:size(d,1)-1,n) = interp1([1:1:size(d0,1)],d0(:,n),[1.5:1:size(d0,1)])';
    p = polyfit([1:3]',d(end-3:end-1,n),2);
    d(end,n) = polyval(p,4);
    p = polyfit([1:3]',d(2:4,n),2);
    d(1,n) = polyval(p,0);
end

