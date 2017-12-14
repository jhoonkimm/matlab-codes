function [output] = caloricEquivalent(oxygen,RER,varargin)
% [output] = caloricEquivalent(oxygen,rer,varargin)
% Used to find caloric equivalents of VO2 data from standard metabolic carts and uses RER to calculate caloric equivalence.
%
% Flags;
% Non-existent currently
%
% Created by Jay Kim (11/11/2017)

%% Default input arguments


%% Optional input arguments

% for i = 1:2:length(varargin)
%   switch varargin{i}
%     case 'eg'
%     eg = varargin{i + 1};
%   end
% end

%% Conversion table

% 1st column is RQ and 2nd column is kJ/L O2 equivalent
enConv = [0.7036	20.287
          0.705     20.291
          0.71      20.316
          0.715     20.341
          0.72      20.366
          0.725     20.387
          0.73      20.412
          0.735     20.437
          0.74      20.463
          0.745     20.488
          0.75      20.509
          0.755     20.534
          0.76      20.559
          0.765     20.584
          0.77      20.605
          0.775     20.63
          0.78      20.655
          0.785     20.68
          0.79      20.705
          0.795     20.726
          0.8       20.751
          0.805     20.776
          0.81      20.801
          0.815     20.826
          0.82      20.847
          0.825     20.872
          0.83      20.897
          0.835     20.923
          0.84      20.943
          0.845     20.969
          0.85      20.994
          0.855     21.019
          0.86      21.044
          0.865     21.065
          0.87      21.09
          0.875     21.115
          0.88      21.14
          0.885     21.161
          0.89      21.186
          0.895     21.211
          0.9       21.236
          0.905     21.261
          0.91      21.282
          0.915     21.307
          0.92      21.332
          0.925     21.357
          0.93      21.378
          0.935     21.403
          0.94      21.429
          0.945     21.454
          0.95      21.479
          0.955     21.5
          0.96      21.525
          0.965     21.55
          0.97      21.575
          0.975     21.596
          0.98      21.621
          0.985     21.646
          0.99      21.671
          0.995     21.696
          0.996     21.7];
      
%% Finding caloric equivalence

cal = [];
   for i = 1:length(RER)
          r = (round(RER(i)*100))/100;
          if r < 0.7050 %lower caloric end
              r = 0.7036;
          elseif r > 0.995 %higher caloric end
              r = 0.996;
          end
      [rw,~] = find(enConv==r);
        cal = [cal;oxygen(i)*enConv(rw,2)*(1/60)]; %gives W kg-1
   end
 output = cal;