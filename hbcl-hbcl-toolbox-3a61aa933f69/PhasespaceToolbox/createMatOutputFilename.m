function [outputFilename] = createMatOutputFilename(inputFilename)
% dotIndex = findstr(inputFilename, '.');
% outputFilename = [inputFilename(1:(dotIndex(end) - 1)) '.mat'];
outputFilename = [inputFilename '.mat'];
end

