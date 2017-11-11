function [varargout] = anyfit(xin,yin,varargin)

% ANYFIT linear fit with offset for multiple subjects
%
% anyfit(xdata,ydata) tries one particular fit for the function
% ydata = c*f(xdata) + d(i),  (default f(xdata) = [xdata])
% *allowing a different offset for each subject's data.
% 
% NOTE that the option "function" overrides the default function with the
% one supplied (see below).  
% Then the fit is:  ydata = c(1)*fn(1) + c(2)*fn(2) + ... + d(i)
% 
% Inputs:  xdata and ydata can be either cell arrays or standard arrays:
% (1) xdata is a cell array, {# of subjects} by (# of trials)
%     ydata is a cell array, {# of subjects} by (# of trials)
% (2) xdata and ydata are same-size standard arrays: 
%       # subjects-by-# trials
%     If standard arrays are used, missing data can be indicated by 
%     NAN values in the "ydata" array; these will be ignored for the fit.  
% 
% Optional:
%     The following options can be appended with commands of the
%     form anyfit(xdata,ydata,'type', 'indiv')
%     'type' = 1 or 'indiv' for individualized fits
%            = 2 or 'single' for single fit (default)
%                individualized fits are where the d(i) is different for each
%                subject.  single fit is where there is only one d(i)
%     'offset' = -1 if true offset is to be restored in plot (default)
%              = y-intercept if needed
%              = 'none' if the fit is to be performed without any
%                y-intercept fit at all
%     'symbol' = the symbol to plot the data points with
%     'axes' = handle to the axes in which to draw the plot
%     'function' = anonymous function and function definition with which to
%               build elements of the prediction matrix from evaluation of the 
%               independent variable.  For example, 
% 
%                   any(xdata,ydata,'function',@(x)[sin(x), cos(x)]) 
%               
%               will return coefficients of sin(x), cos(x), and y-offset (1 by default)
%                   (i.e., y = a*sin(x) + b*cos(x) + c)
%               that best fit the y vs x data.  Other options should still work 
%               as specified above.  Specifically, for a "no offset" fit 
%               ( y = a*sin(x) + b*cos(x) ), use options {'offset','none'}.
%     'colors' = matrix of colors representing the color order to use for
%               different subjects' data
%     'alpha' = Alpha value for the 100*(1-Alpha)% confidence interval
%
% Outputs:
%       coef    = [c(:); d(:)]
%       r2      is the R^2 value for the fit
%       coefint is the 95% confidence interval for each coef
%       stats   is the regression R^2, F-value, and p-value
%       xd      is the x data with missing data points removed
%       yd      is the y data with missing data points removed and offsets adjusted
%       offset  is the overall offset used for the curve fit.  Equals the
%               mean of the offset coefficients for each subject for 'indiv' type
%               fits, or the global offset for 'single' type fits.  
%       hfig    is the figure handle for the resulting plot
%       hax     is the axes handle for the resulting plot
%       hdata   is the handle array to the data points
%       hline   is the handle to the fit line
% 
%       If only one output is requested, the function returns a structure
%       with all the above fields.
%
% Also, the mean data plus error bars are plotted, along with the
% fit.  All data are adjusted so that the fit has a y-intercept of
% zero or offset if that input is provided.  An offset of -1 means
% that even for an individualized fit, the true offset in the data
% is restored


% Art Kuo 12/2000, latest version has optional arguments
% remade for general use by Peter Adamczyk, 11/2005 (formerly "artfit.m")


% handle inputs for cases of cell array or not.
xdata = {};
ydata = {};
if iscell(xin)
    if iscell(yin)
        xdata = xin;
        ydata = yin;
%         xdata = reshape(xin,1,[]) ;
%         ydata = reshape(yin,1,[]) ;
    else
        error('data type mismatch between x and y: x is a cell array and y is not.')
        return
    end
else
    if iscell(yin)
        error('data type mismatch between x and y: y is a cell array and x is not.')
        return
    end
    % if we get here, then neither xin nor yin is a cell array.  Handle them as
    % arrays.  
    xdata = num2cell(xin,2);
    ydata = num2cell(yin,2);
%     for i = 1:size(yin,1)
%         % a value of 0 or NaN in the Y-data indicates nonexistent data point
%         goods = find( (yin(i,1:end) ~= 0) & ~isnan(yin(i,1:end)) );
%         xdata{i} = xin(i,goods);
%         ydata{i} = yin(i,goods);
%     end
end
xdata = reshape(xdata,1,[]) ;
ydata = reshape(ydata,1,[]) ;

% Each ROW has one subject's data, and NaN values indicate missing data.
%% removed by Peter 2007/05/30 because you can leave the NaN's in and it
%% causes no problems*.   
%% * Thank goodness the calculation of R^2 can use Sum Squared Error, and
%% does not require Mean Square Error! 
% for i = 1:length(xdata)
%     goods = find( ~isnan(xdata{i}) & ~isnan(ydata{i}) );
%     xdata{i} = xdata{i}(goods);
%     ydata{i} = ydata{i}(goods);
% end

        
% data is now a cell array of # of subjects  
nsubs = length(xdata);

if length(xdata) ~= length(ydata)
  error('xdata and ydata do not match')
end

% The A matrix depends on the type of fit.  Its general form is
% A = [f(x) 1].  The 1 matrix depends on the type of fit.  
% If type = 1, 1 is a matrix that is ntrials*nsubs by nsubs,
% and it contains columns of ntrials 1's, each staggered from
% another.  If type = 2, 1 is just a long vector of 1's.

varargin = reshape(varargin,2,[])'; 
i = strmatch('type',{varargin{:,1}});
if ~isempty(i) & ...
  (strncmpi('indiv',varargin{i,2},5) | varargin{i,2}==1),
  ncols = nsubs; type = 1;
elseif isempty(i) | strmatch('single',varargin{i,2})
  ncols = 1; type = 2;
end

% even if A _should_ include a "1" matrix, it can be forced not to by the
% 'offset','none' option.  Enforced by setting "ncols" to 0 (width of the
% "1" portion of the A matrix.
i = strmatch('offset',{varargin{:,1}});
if ~isempty(i)
  offset = varargin{i,2};
  if strmatch('none',offset)  % fit everything without an offset 
    offset = 0;
    ncols = 0;
  end
else
  offset = -1;
end

i = strmatch('axes',{varargin{:,1}});
if ~isempty(i)
  axhandle = varargin{i,2};  
else
  axhandle = [];
end

i = strmatch('symbol',{varargin{:,1}});
if ~isempty(i)
  s = varargin{i,2};
else
  s = '.';
end

i = strmatch('alpha',{varargin{:,1}});
if ~isempty(i)
  alpha = varargin{i,2};
else
  alpha = 0.05 ;
end

% Arrange X-data in nice, long columns
xvalues = cat(2,xdata{:})'; 
xcols = size(xvalues,2);

i = strmatch('function',{varargin{:,1}});
if ~isempty(i)
  funchandle = varargin{i,2};  % the "@(x) [sin(x) cos(x)]" syntax _may_ need to be revised for Matlab R15.  Maybe.  If so, move to CELL array.  
  nfuncs = length(funchandle(ones(1,xcols)));
else
  funchandle = @(x) [x];
  nfuncs = length(funchandle(ones(1,xcols)));
end

% it is assumed that each subject has their own set of x-values
% xvalues = cat(2,xdata{:})'; 
% xcols = size(xvalues,2);
A = [feval(funchandle,xvalues) zeros(length(xvalues), ncols)];

% Now make the columns of 1's for the "offset" computation:
ntrialtot = 0;
if type == 1 
  for i=1:ncols
    ntrials = length(ydata{i});
    A(ntrialtot+(1:ntrials),nfuncs+i) = 1;
    ntrialtot = ntrialtot + ntrials;
  end
elseif ncols == 1
  A(:,nfuncs+1) = 1;
end

y = cat(2,ydata{:})';

% get rid of any row containing a NaN  
badrows = union(find(any(isnan(A),2)),find(any(isnan(y),2)));
A(badrows,:)=[]; y(badrows) =[]; xvalues(badrows)=[];
badcols = find(~any(A,1)); % columns with no Offset column remaining (e.g., whole subject removed or empty) are no good for the offset coefficient averaging

coef = A\y ;
% xstart = [coef];%; wn]; 

% R^2 calculation
sst=sum((y-mean(y)).^2);
sse=sum((y-A*coef).^2);
r2=1-sse/sst;
% fprintf(1,'r2 = %g\n', r2);

coef(badcols) = NaN;%% Offset columns that were empty need to have a NaN in the coefficient, otherwise they will mess up the Offset Mean

[b,bint,r,rint,stats] = regress(y,A,alpha);
coefint = bint(:,2)-b;  % "bint" is actually the 95% confidence limit. Subtract the nominal value "b" to get a +/- margin.
coefint(badcols) = NaN;
%clf; plot(r); pause


% produce a smaller mesh of x to make a smoother line
xmin = min(xvalues,[],1) ;
xmax = max(xvalues,[],1) ;
xs(:,1) = linspace(xmin(1),xmax(1),30);
% To make this work for multi-column inputs, fix this next part.
% 
% % [xvb,sortind] = sort(xvalues,1,'ascend') ;
% % xvb = xvalues(sortind,:);
% % for i = 2:xcols
% % %     xs(:,i);
% %     q = smooth(xvb(:,1),xvb(:,i),10,'rlowess'); 
% %     xs(:,i) = interp1(xvb(:,1),q,xs(:,1));
% % %     size(xs)
% % %     size(q)
% % end
% % 

if offset == -1
  offsetval = nanmean(coef(nfuncs+1:end));
else 
    offsetval = 0;
end


for i=1:nsubs
  if type == 1 & offset == -1 % individualized fit
    fixy = coef(i+nfuncs);  % subtract the individual offset to each subject's data
  elseif ncols == 1 & offset == -1 % single fit, non-zero offset
    fixy = coef(nfuncs+1);    % subtract the single offset to each subject's data
  else              % single fit, zero offset
    fixy = 0;          % add nothing to each subject's data
  end
  ydataoff{i} = ydata{i} - fixy + offsetval;
  syms{i} = s; 
end


% % % % % Plotting of the graph
if isempty(axhandle)
    hfig = figure; clf;
    axhandle = axes;
else
    axes(axhandle);
    hfig = gcf;
end

% colors specified?
i = strmatch('colors',{varargin{:,1}});
if ~isempty(i)
    colors = varargin{i,2}; 
    if iscell(colors)
        colors = cat(1, colors{:});
    end
  set(axhandle, 'ColorOrder',colors);  
else
  %do nothing
end

% arcfootmap = [get(gcf,'DefaultAxesColorOrder'); 1 .5 0; .5 0 .5; .25 .85 .25; 0 0 0]; arcfootmap(5,:) = [1 0 1];
% set(gcf,'DefaultAxesColorOrder',arcfootmap);
hline = plot(axhandle,xs,feval(funchandle,xs)*coef(1:(nfuncs)) + offsetval,'b-','LineWidth',3); hold on % plot the fit

crap = [xdata; ydataoff; syms]; % plot the data points
hdata = plot(axhandle, crap{:});
set(hdata,'MarkerSize',10)
v = axis; 

axis([0 v(2) 0 v(4)*1.2]);

funcstring = func2str(funchandle) ;
text(.1*v(2),1.1*v(4), {[funcstring(5:end) ' * [' num2str(coef(1:nfuncs)', '%10.5f' ) ']'' + ' num2str(offsetval,'%10.5f' ) ], ['CI +/- ' num2str([coefint(1:nfuncs)', nanmean(coefint(nfuncs+1:end))], '%10.5f')], ['R^2 value = ' num2str(r2,'%5.4f') ]})

%%%%%%  Print to Screen
    fprintf(1,'r2 for given fit = %g\n', r2);
    fprintf(1,'For coefficients \n')
    fprintf(1,'         %g\n', coef(1:nfuncs), offsetval );
%%%%%%%%%


%errorbar(xdata,ydata
%errorbar(mean(x,2),mean(data-repmat(coef(p+1:end)',ntrials,nsubs/ncols),2)+offset,...
%  std(data-repmat(coef(p+1:end)',ntrials,nsubs/ncols),[],2),'r.')

hold off

yd = ydataoff;
xd = xdata;


% distribute outputs
[var{1:8}] = deal(coef,r2,coefint,stats,xd,yd,offset,hfig);
% if nargout == 1
    vars.coef = coef;
    vars.r2 = r2;
    vars.coefint = coefint;
    vars.stats = stats;
    vars.xd = xd;
    vars.yd = yd;
    vars.offset = offsetval;
    vars.hfig = hfig;
    vars.hax = axhandle;
    vars.hdata = hdata;
    vars.hline = hline;
if nargout == 1
    varargout{1} = vars;
else
    for i = 1:nargout
            varargout{i} = var{i};
    end
end
    
set(axhandle, 'UserData', vars);
drawnow


%% The following is a Fix for Individual Subjects Fit, to compute the R^2
%% value AFTER shifting all the data by the individual offsets
% if 1 & (type==1) & (offset == -1)
%     coeffix = coef; coeffix(1) = 0;
%     yfix = y - A*coeffix + mean(coeffix(2:end));
%     coeffix = A\yfix;
%     sstfix = sum((yfix-mean(yfix)).^2);
%     ssefix = sum((yfix-A*coeffix).^2);
%     r2fix=1-ssefix/sstfix;
%     [bfix,bintfix,rfix,rintfix,statsfix] = regress(yfix,A,alpha);
%     coefintfix = bintfix(:,2)-bfix;  % "bint" is actually the 95% confidence limit. Subtract the nominal value "b" to get a +/- margin.
%     coefintfix(badcols) = NaN;
% end
% 
% % Label the graph with "fixed" parameter values and R^2, after accounting
% % for the effect of individual offsets. 
% text(.1*v(2),1.1*v(4), {[funcstring(5:end) ' * [' num2str(coeffix(1:nfuncs)', '%10.5f' ) ']'' + ' num2str(offsetval,'%10.5f' ) ], ['CI +/- ' num2str([coefintfix(1:nfuncs)', nanmean(coefintfix(nfuncs+1:end))], '%10.5f')], ['R^2 value = ' num2str(r2fix,'%5.4f') ]})
