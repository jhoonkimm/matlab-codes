function [pvalues, results] = hbclstats(data,ctrl,alpha);
% performs stats on a repeated-measure table of data, comparing all
% columns against each other. Each row stands for one subject.
% This assumes you have already performed an appropriate ANOVA
% with repeated measures. Assuming it shows significance, then
% use this to perform a set of paired t-tests under the 
% Holm-Sidak step-down procedure.
%
% hbclstats(data,ctrl) where ctrl = N uses the Nth column as a control,
%   and performs only comparisons against control. When ctrl = 0, all
%   pair-wise comparisons are made.
% hbclstats(data,ctrl,alpha) uses the designated value for alpha.
%
% [pvalues, results] = hbclstats(...) returns the p-values in a list
%   and the list of results in a step-down table. The pvalues list
%   is in the order of comparisons. For all comparisons, that order
%   is 1-2, 1-3, ..., 1-n, 2-3, 2-4, ...  For comparisons against a
%   control, the order is ctrl-1, ctrl-2, ctrl-3, ... (note that the
%   list is n-1 elements long).
%   The results table has columns v1, v2, pvalue, alpha, reject
%   where v1 and v2 are the values being compared, alpha is the 
%   corrected alpha for step-down comparison, and reject is 1 if
%   the null hypothesis can be rejected.
%
% Primer of Biostatistics
% By Stanton A. Glantz
% Published by McGraw-Hill Professional, 2005

% Art Kuo 3/2008  

if nargin == 1,
  ctrl = 0;
end
if nargin <= 2,
  alpha = 0.05;
end

if isempty(ctrl), 
  ctrl = 0;
end

[m,n] = size(data);

pvalues = ones(n*(n-1)/2, 1)*NaN;
comparolist = zeros(n*(n-1)/2, 2); 

k = 0;

if ctrl == 0, % perform t-tests on everything, put them in a list of pvalues
  for i = 1:n-1
    for j = i+1:n
      k = k + 1;
      [h, pvalues(k)] = ttest(data(:,i), data(:,j));
      comparolist(k, 1:2) = [i j];
    end
  end
else % perform tests against a control
  pvalues = ones(n-1, 1)*NaN; comparolist = zeros(n-1, 2);
  for j = setdiff(1:n, ctrl)
    k = k + 1;
    [h, pvalues(k)] = ttest(data(:,ctrl), data(:,j));
    comparolist(k, 1:2) = [ctrl j];
  end
end

[pvaluessorted, sortlist] = sort(pvalues);
stepdownlist = 1 - (1-alpha).^(1./(k:-1:1)');

results = [comparolist(sortlist,:) pvaluessorted stepdownlist pvaluessorted<stepdownlist];

fprintf(1, 'Holm-Sidak procedure results\n\n');
fprintf(1, 'test\tp-values\talpha\t reject\n');
fprintf(1, '% d-%d\t %6.2g \t%6.2g\t  %d\n', results');
fprintf(1, '\n');
