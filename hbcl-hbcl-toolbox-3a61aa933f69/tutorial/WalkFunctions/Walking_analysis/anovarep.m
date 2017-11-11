function pvalues = anovarep(varargin)
% ANOVAREP - Analysis of variance for repeated measures. 
% This function executes the analysis of variance when subjects underwent
% several treatments. This function is similar to ANOVA2 Matlab Function,
% but there are three differences:
% 1) the output of ANOVA table;
% 2) the graphical plot of anova table;
% 3) if p-value<alpha this function executes the Holm-Sidak test for multiple
% comparison test to highlight differences between treatments.
% 
% Syntax: 	ANOVAREP(X,ALPHA [, DISPTABLE])
%      
%     Inputs:
%           X - data matrix. 
%           ALPHA - significance level (default = 0.05).
%           DISPTABLE - set to 1 to produce graphical tables
%     Outputs:
%           - Anova table.
%           - Graphical plot of anova table.
%           - Holm-Sidak test (eventually)
%
%     pvalues = anovarep(x, alpha) returns a structure with
%           elements anova and table, containing results of
%           anova (between subjects and between columns) and
%           the table of Holm-Sidak values
%
%      Example: demoanovarep
%
%           Created by Giuseppe Cardillo
%           giuseppe.cardillo-edta@poste.it
%
% To cite this file, this would be an appropriate format:
% Cardillo G. (2008) Anovarep: compute the Anova for repeated measures and
% Holm-Sidak test for multiple comparisons if Anova is positive. 
% http://www.mathworks.com/matlabcentral/fileexchange/18746

% This version has an option to display the table or not (comment added by
% Art)

%Input error handling
args=cell(varargin);
nu=numel(args);
if isempty(args)
    error('stats:anovarep:TooFewInputs','At least one input is required.');
elseif nu>3 
    error('stats:anovarep:TooMuchInputs','Max three inputs are required.');
else
% following line modified by Art to include third argument
    default.values = {[],0.05, 1}; %default value of alpha
    default.values(1:nu) = args;
% and another modification by Art to include disptable
    [x alpha disptable] = deal(default.values{:});
    if isvector(x)
        error('Warning: x must be a matrix, not a vector.');
    end
    if (any(isnan(x(:)))) %check the matrix
        error('stats:anovarep:DataNotBalanced',...
            'NaN values in input not allowed.  Use anovan instead.');
    end
    if ~all(isnumeric(x(:))) || any(isinf(x(:))) || isempty(x)
        error('The x values must be numeric and finite')
    end
    if nu==2 %if alpha was given
        if ~isscalar(alpha) || ~isnumeric(alpha) || ~isfinite(alpha) || isempty(alpha)
            error('Warning: ALPHA value must be a numeric and finite scalar');
        end
        if alpha <= 0 || alpha >= 1 %check if alpha is between 0 and 1
            error('Warning: ALPHA must be comprised between 0 and 1.')
        end        
    end
    if (isnan(alpha) || alpha <= 0 || alpha >= 1)%check if alpha is between 0 and 1
        error('Warning: ALPHA must be comprised between 0 and 1.')
    end        
end

clear args default nu

[n,m]=size(x); %n=subjects; m=treatments;
T=mean(x); %Mean of treatment
S=mean(x,2); %Mean of subject
Xbar=sum(x(:))/(m*n); %General mean

SStra=m*sum((S-Xbar).^2); %Variability among subjects
GLtra=n-1; %degrees of freedom
MStra=SStra/GLtra; % variance of population estimed from subjects

SSentro=sum(sum((x-repmat(S,1,m)).^2,2)); %Variability within subjects
GLentro=n*(m-1); %degrees of freedom

SStot=SStra+SSentro; %Total Variability
GLtot=GLtra+GLentro; %Degrees of freedom

%The variability within subjects can be decomposed in:

SSt=n*sum((T-Xbar).^2); %Variability among treatments
GLt=m-1; %degrees of freedom
MSt=SSt/GLt; % variance of population estimed from treatments

SSr=SSentro-SSt; %Residual variability
GLr=(n-1)*(m-1); %degrees of freedom
MSr=SSr/GLr; % variance of population estimed from residuals

F=[MStra MSt]./MSr;
panova=1-fcdf(F,[GLtra GLt],GLr); %p-value

%Formatting for ANOVA Table printout.
Table{6,6} = [];   
Table(1,1:6) = {'Variability' 'SS' 'df' 'MS' 'F' 'p-value'};
Table(2:6,1)={'Total';'Among subjects';'Within subjects';'Among treatments';'Residual'};
Table(2:6,2)={SStot; SStra; SSentro; SSt; SSr};
Table(2:6,3)={GLtot; GLtra; GLentro; GLt; GLr};
Table([3 5 6],4)={MStra; MSt; MSr};
Table([3 6],5)={F(1); F(2)};
Table([3 6],6)={panova(1); panova(2)};

if panova(2)<alpha
    footer='Almost one of the treatments is the cause of variation';
else
    footer='Variation is due to chance';
end

if disptable % addition by Art
  H=figure;%('Position',[20 803 560 147]);
  statdisptable(Table, 'ANOVA FOR REPEATED MEASURES', 'ANOVA Table', footer,[-1 -1 0 -1 2 4],H);

  %Display bar a bar graph of anova table
  H=figure;
  %set(H,'Position',[600 502 560 420]);
  y=[0 0 1 1];
  hold on
  h1=fill([0 SStot SStot 0],y,'c');
  y=y+1;
  h2=fill([0 SStra SStra 0],y,'y');
  h3=fill([0 SSentro SSentro 0]+SStra,y,'r');
  y=y+1;
  h4=fill([0 SSt SSt 0]+SStra,y,'g');
  h5=fill([0 SSr SSr 0]+SStra+SSt,y,'b');
  hold off
  legend([h1 h2 h3 h4 h5],'Total','Among subjects','Within subjects','Among treatments','Residual','Location','NorthWest')
  title('Sources of Variability')
  clear S Xbar SStra GLtra SSentro GLentro SStot GLtot SSt GLt MSt SSr F
  clear y h1 h2 h3 h4 h5 Table
else
  clear Table
end

%If Panova<alpha execute the Holm-Sidak test to close off the differences
%in anova
if panova(2)<alpha
    a=m-1; %rows of probability matrix
    c=0.5*m*(m-1); %max number of comparisons
    count=0; %counter
    p=ones(1,c); %preallocation of p-value vector
    pb{c,2} = []; %preallocation of p-value matrix
    denom=realsqrt(2*MSr/n); %the denominator is the same for each comparison
    for I=1:a
        for J=I+1:m
            count=count+1;
            t=abs(diff(T([I J])))/denom; %t-value
            p(count)=(1-tcdf(t,GLr))*2; %2-tailed p-value vector
            pb(count,:)={strcat(int2str(I),'-',int2str(J));num2str(p(count))}; %Matrix of the p-values
        end
    end
    clear a m t count I J %clear unnecessary variables
    [p,I]=sort(p); %sorted p-values
    pb=pb(I,:);
    J=1:c; %How many corrected alpha value?
    alphacorr=1-((1-alpha).^(1./(c-J+1))); %Sidak alpha corrected values
    %Compare the p-values with alpha corrected values. 
    %If p<a reject Ho; else don't reject Ho: no more comparisons are required.
    Table{c+1,4} = [];   %Formatting for Holm-Sidak Table printout.
    Table(1,1:4) = {'Test' 'p-value' 'alpha' 'Comment'};
    comp=1; %compare checker
    row=2;
    Table(2:end,1:2)=pb;
    for J=1:c
        if comp %Comparison is required
            Table(row,3)={num2str(alphacorr(J))};
            if p(J)<alphacorr(J)
                Table(row,4)={'Reject H0'};
            else
                Table(row,4)={'Fail to reject H0'};
                comp=0; %no more comparison are required
            end
        else %comparison is unnecessary
            Table(row,3:4)={'No comparison made';'H0 is accepted'};
        end
        row=row+1;
    end
    if disptable
      H=figure;%('Position',[20 320 560 420]);
      statdisptable(Table,'','Holm-Sidak multiple t-test','',[-1 -1 0 -1 2 4],H);
    end
else
  p = NaN;
end

pvalues.anova = panova;

pvalues.table = p;
