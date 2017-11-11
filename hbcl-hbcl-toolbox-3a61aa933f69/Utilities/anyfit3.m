function [varargout] = anyfit3(xin,yin,zin,varargin)

% ANYFIT3 3D linear fit with offset for multiple subjects
%
% anyfit3(xdata,ydata,zdata) tries one particular fit for the function
% zdata = c*f(xdata) + d*g(ydata) + e(i),  
% (default f(xdata) = [xdata] and f(ydata) = ydata)
% *allowing a different offset for each subject's data.
% 
% NOTE that the option "function" overrides the default function with the
% one supplied (see below).  
% Then the fit is:  zdata = c(1)*fn(1) + c(2)*fn(2) + ... + d(i)
% 
% Inputs:  xdata, ydata and zdata can be cell arrays or standard arrays:
% (1) xdata is a cell array, {# of subjects} by (# of trials)
%     ydata is a cell array, {# of subjects} by (# of trials)
%     zdata is a cell array, {# of subjects} by (# of trials)
% (2) xdata and ydata and zdata are same-size standard arrays: 
%       # subjects-by-# trials
%     If standard arrays are used, missing data can be indicated by 
%     NAN values in the "ydata" array; these will be ignored for the fit.  
% 
% Optional:
%     The following options can be appended with commands of the
%     form anyfit(xdata,ydata,zdata,'type', 'indiv')
%     'type' = 1 or 'indiv' for individualized fits
%            = 2 or 'single' for single fit (default)
%                individualized fits are where the d(i) is different for each
%                subject.  single fit is where there is only one d(i)
%     'offset' = -1 if true offset is to be restored in plot (default)
%              = z-intercept if needed
%              = 'none' if the fit is to be performed without any
%                z-intercept fit at all
%     'symbol' = the symbol to plot the data points with
%     'axes' = handle to the axes in which to draw the plot
%     'function' = anonymous function and function definition with which to
%               build elements of the prediction matrix from evaluation of the 
%               independent variable.  For example, 
% 
%                   any(xdata,ydata,zdata,'function',@(x,y)[sin(x), cos(y)]) 
%               
%               will return coefficients of sin(x), cos(y), and z-offset (1 by default)
%                   (i.e., z = a*sin(x) + b*cos(y) + c)
%               that best fit the z, x and y data.  Other options should still work 
%               as specified above.  Specifically, for a "no offset" fit 
%               ( z = a*sin(x) + b*cos(y) ), use options {'offset','none'}.
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
% fit.  All data are adjusted so that the fit has a z-intercept of
% zero or offset if that input is provided.  An offset of -1 means
% that even for an individualized fit, the true offset in the data
% is restored


% Art Kuo 12/2000, latest version has optional arguments
% remade for general use by Peter Adamczyk, 11/2005 (formerly "artfit.m")
% modified from "anyfit.m" by Peter Adamczyk 01/2008 for 3D fits


% handle inputs for cases of cell array or not.
xdata = {};
ydata = {};
zdata = {};
if iscell(zin)
    if iscell(yin)&iscell(xin)
        xdata = xin;
        ydata = yin;
        zdata = zin;
%         xdata = reshape(xin,1,[]) ;
%         ydata = reshape(yin,1,[]) ;
    else
        error('data type mismatch between x, y, z: z is a cell array and x or y is not.')
        return
    end
else % z is NOT a cell array...
    if iscell(yin)|iscell(xin)
        % but x or y is
        error('data type mismatch between x, y, z: x or y is a cell array and z is not.')
        return
    end
    % if we get here, then neither nothing is a cell array.  Handle them as
    % arrays.  
    xdata = num2cell(xin,2);
    ydata = num2cell(yin,2);
    zdata = num2cell(zin,2);
%     for i = 1:size(yin,1)
%         % a value of 0 or NaN in the Y-data indicates nonexistent data point
%         goods = find( (yin(i,1:end) ~= 0) & ~isnan(yin(i,1:end)) );
%         xdata{i} = xin(i,goods);
%         ydata{i} = yin(i,goods);
%     end
end
xdata = reshape(xdata,1,[]) ;
ydata = reshape(ydata,1,[]) ;
zdata = reshape(zdata,1,[]) ;

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

if (length(zdata) ~= length(ydata)) | (length(zdata) ~= length(xdata))
  error('zdata and ydata or xdata do not match')
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

% Arrange Y-data in nice, long columns
yvalues = cat(2,ydata{:})'; 
ycols = size(yvalues,2);

% is a new "function" supplied?
i = strmatch('function',{varargin{:,1}});
if ~isempty(i)
  funchandle = varargin{i,2};  % the "@(x) [sin(x) cos(x)]" syntax _may_ need to be revised for Matlab R15.  Maybe.  If so, move to CELL array.  
  nfuncs = length(funchandle(ones(1,xcols),ones(1,ycols)));
else
  funchandle = @(x,y) [x y];
  nfuncs = length(funchandle(ones(1,xcols),ones(1,ycols)));
end

% it is assumed that each subject has their own set of x- and y-values
% xvalues = cat(2,xdata{:})'; 
% yvalues = cat(2,ydata{:})'; 
% xcols = size(xvalues,2);
% ycols = size(yvalues,2);
A = [feval(funchandle,xvalues,yvalues) zeros(length(xvalues), ncols)];

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

z = cat(2,zdata{:})';

% get rid of any row containing a NaN  
badrows = union(find(any(isnan(A),2)),find(any(isnan(z),2)));
A(badrows,:)=[]; z(badrows) =[]; xvalues(badrows)=[]; yvalues(badrows)=[];
badcols = find(~any(A,1)); % columns with no Offset column remaining (e.g., whole subject removed or empty) are no good for the offset coefficient averaging

coef = A\z ;
% xstart = [coef];%; wn]; 

sst=nansum((z-mean(z)).^2);
sse=nansum((z-A*coef).^2);
r2=1-sse/sst;
% fprintf(1,'r2 = %g\n', r2);

coef(badcols) = NaN;%% Offset columns that were empty need to have a NaN in the coefficient, otherwise they will mess up the Offset Mean


[b,bint,r,rint,stats] = regress(z,A,alpha);
coefint = bint(:,2)-b;  % "bint" is actually the 95% confidence limit. Subtract the nominal value "b" to get a +/- margin.
coefint(badcols) = NaN;
%clf; plot(r); pause


% produce a smaller mesh of x and y to make a smoother line
xmin = min(xvalues,[],1) ;
xmax = max(xvalues,[],1) ;
xs(:,1) = linspace(xmin(1),xmax(1),30);
ymin = min(yvalues,[],1) ;
ymax = max(yvalues,[],1) ;
ys(1,:) = linspace(ymin(1),ymax(1),30);

xs = repmat(xs,1,size(ys,2));
ys = repmat(ys,size(xs,1),1);

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
    fixz = coef(i+nfuncs);  % subtract the individual offset to each subject's data
  elseif ncols == 1 & offset == -1 % single fit, non-zero offset
    fixz = coef(nfuncs+1);    % subtract the single offset to each subject's data
  else              % single fit, zero offset
    fixz = 0;          % add nothing to each subject's data
  end
  zdataoff{i} = zdata{i} - fixz + offsetval;
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
zmesh = reshape(feval(funchandle,xs(:),ys(:))*coef(1:(nfuncs)),size(xs,1),size(xs,2)) + offsetval ;
hline = surf(axhandle,xs,ys,zmesh,zmesh,'FaceAlpha',0.4,'EdgeAlpha',0,'FaceColor','interp','EdgeColor','interp'); hold on % plot the fit
% hline = mesh(axhandle,xs,ys,zmesh,zmesh); hold on % plot the fit

crap = [xdata; ydata; zdataoff; syms]; % plot the data points
hdata = plot3(axhandle, crap{:});
set(hdata,'MarkerSize',10)
v = axis; 

axis([v(1:5) v(6)*1.2]);

funcstring = func2str(funchandle) ;
text(v(1)+.1*(v(2)-v(1)),v(3)+.1*(v(4)-v(3)),v(5)+1.1*(v(6)-v(5)),{[funcstring ' * [' num2str(coef(1:nfuncs)', '%10.5f' ) ']'' + ' num2str(offsetval,'%10.5f' ) ], ['CI +/- ' num2str([coefint(1:nfuncs)', nanmean(coefint(nfuncs+1:end))], '%10.5f')], ['R^2 value = ' num2str(r2,'%5.4f') ]})

%%%%%%  Print to Screen
    fprintf(1,'r2 for given fit = %g\n', r2);
    fprintf(1,'For coefficients \n')
    fprintf(1,'         %g\n', coef(1:nfuncs), offsetval );
%%%%%%%%%


%errorbar(xdata,ydata
%errorbar(mean(x,2),mean(data-repmat(coef(p+1:end)',ntrials,nsubs/ncols),2)+offset,...
%  std(data-repmat(coef(p+1:end)',ntrials,nsubs/ncols),[],2),'r.')

hold off

zd = zdataoff;
xd = xdata;
yd = ydata;


% distribute outputs
[var{1:9}] = deal(coef,r2,coefint,stats,xd,yd,zd,offset,hfig);
% if nargout == 1
    vars.coef = coef;
    vars.r2 = r2;
    vars.coefint = coefint;
    vars.stats = stats;
    vars.xd = xd;
    vars.yd = yd;
    vars.zd = zd;
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