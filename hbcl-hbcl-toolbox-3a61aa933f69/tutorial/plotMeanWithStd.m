function h=plotMeanWithStd(x,y,lin,thickness,trans,transcol)
% h=plotMeanWithStd(x,y,lin,thickness,trans,transcol)
% plots y and +/- standard deviation of y against x
%   x should be a row vector
%   y should be a matrix of row vectors
%   lin is the line type for the plot of y in standard plot() format
%   thickness is the thickness value for that line (ex. 0.2)
%   trans is the transparency amount (0 to 1)
%   transcol is the color of the band in [r g b] where each is a value from
%       0 to 1
% 
poly=mean(y)-std(y);
hfill=fill([x x(end:-1:1)],[mean(y)+std(y) poly(end:-1:1)],transcol);
set(hfill,'EdgeColor',transcol)
% set(hfill,'FaceColor',transcol)
set(hfill,'FaceAlpha',trans)
hold on
h=plot(gca,x,mean(y),lin,'LineWidth',thickness);
% h=plot(gca,x,mean(y),line,x,mean(y)+std(y),[line(1) ':'],x,mean(y)-std(y),[line(1) ':'],'LineWidth',thickness);
% h=h(1);
end