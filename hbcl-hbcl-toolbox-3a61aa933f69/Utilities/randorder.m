function neworder = randorder(n)

% RANDORDER randomly reorder list
% Argument N
%  -if N is integer, this reorders list from 1:n
%  -if N is an array, this reorders array values
%
% Example:
%  randorder(5)
%  randorder([10 20 30 40])

% Karl Zelik 4/13/10

% Obsolete way
% % % function neworder = randorder(listn)
% % % RANDORDER randomly reorder a list
% % %  Input array. Outputs array values in random order
% % 
% % % Created by Steve Collins
% % % Updated/documented by Karl Zelik 11/19/09
% % 
% % if length(listn) == 1
% %     neworder = listn;
% % else
% %     rand('state',sum(100*clock)); 
% %     i = round(rand*length(listn) + 0.5);
% %     neworder = [listn(i) randorder([listn(1:i-1) listn(i+1:end)])];
% % end

if length(n) == 1
    A = rand(1,n);
    for i = 1:n
        A(find(A==min(A))) = i;
    end
    
else
    A = rand(1,length(n));
    for i=1:length(n)
        A(find(A==min(A))) = n(i);
    end
end

neworder = A;
    