function [newarray] = nanToZero(oldarray)

% NANTOZERO convert all NaN values in array to zeros
%  
%  Examples:
%   q(1:1000,1)=NaN; q=nanToZero(q);

%  Karl Zelik
%  11/9/09


newarray = oldarray;
i=isnan(newarray);
newarray(i)=0;



% % %% Old (slower) way
% % tic
% % newarray = oldarray;
% % for i=1:size(newarray,1)
% %     for j=1:size(newarray,2)
% %         if isnan(newarray(i,j))
% %             newarray(i,j)=0;
% %         end
% %     end
% % end
% % toc
