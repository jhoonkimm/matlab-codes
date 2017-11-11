function [combinations] = getCombinationsOfListElements(listElements, numberOfElementsInCombination, varargin)
%GETCOMBINATIONSOFLISTELEMENTS Summary of this function goes here
%   Detailed explanation goes here

getSortedUniqueCombinations = 0;

for i = 1:2:length(varargin)
  if (strcmp('getSortedUniqueCombinations', varargin{i}))
    getSortedUniqueCombinations = varargin{i + 1};
  end
end


% vectorsToCombine = cell(numberOfElementsInCombination, 1);
% for i = 1:length(vectorsToCombine)
%   vectorsToCombine{i} = listElements;
% end
% allCombinationsOfIndeces = combvec(vectorsToCombine{:});
% combinations = allCombinationsOfIndeces * NaN;
% indexIntoOutput = 1;
% for i = 1:size(allCombinationsOfIndeces, 2)
%   if (length(unique(allCombinationsOfIndeces(:, i))) == numberOfElementsInCombination)
%     combinations(:, indexIntoOutput) = allCombinationsOfIndeces(:, i);
%     indexIntoOutput = indexIntoOutput + 1;
%   end
% end
% combinations(:, isnan(combinations(1, :))) = [];

if (numberOfElementsInCombination == 1)
  combinations = listElements;
else
  
  combinations = zeros(numberOfElementsInCombination, ...
    nchoosek(length(listElements), numberOfElementsInCombination) * factorial(numberOfElementsInCombination));
  indexIntoCombinations = 0;
  combinations = [];
  for i = 1:length(listElements)
    %     listElements([(1:(i-1)) (i+1):end]);
    %     rest = getCombinationsOfListElements(listElements([(1:(i-1)) (i+1):end]), numberOfElementsInCombination - 1);
    rest = getCombinationsOfListElements(listElements([(1:(i-1)) (i+1):end]), numberOfElementsInCombination - 1);
    theseElements = [repmat(listElements(i), [1, size(rest, 2)]); rest];
    
    %     combinations = [combinations theseElements]
    combinations(:, (1:size(theseElements, 2)) + indexIntoCombinations) = theseElements;
    indexIntoCombinations = indexIntoCombinations + size(theseElements, 2);
  end
  
end



if (getSortedUniqueCombinations)
  combinations = unique(sort(combinations)', 'rows')';
end


end

