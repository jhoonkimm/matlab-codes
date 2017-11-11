function [lmvals,indds]=lminmat(xx,filt,varargin)
%LMINMAT 	function [lmval,indd]=lminmat(x,filt,dim)
%	Find local minima in matrix X, where LMVAL is the output
%	cell array with minima values, INDD is the corresponding indeces 
%	FILT is the number of passes of the small running average filter
%	in order to get rid of small peaks.  Default value FILT =0 (no
%	filtering). FILT in the range from 1 to 3 is usially sufficient to 
%	remove most of a small peaks
%	see also LMIN, MAX, MIN
	
%
%**************************************************|
% 	Serge Koptenko, Guigne International Ltd., |
%	phone (709)895-3819, fax (709)895-3822     |
%--------------06/03/97----------------------------|

%************************************************
% Updated for matrix operations by Peter Adamczyk
% 29 Nov 2007
%************************************************

sz = size(xx);
  if nargin <2
      filt=0; 
      dim = find(sz>1,1,'first') ; % first nonsingleton dimension
  elseif nargin == 2
      dim = find(sz>1,1,'first') ; 
  elseif nargin == 3
      dim = varargin{1};
  else % else too many arguments
      error('Unknown Argument Structure');
  end
  
  xx = num2cell(xx,dim);
  lmvals = cell(size(xx));
  indds = lmvals;
  
  for qq = 1:numel(xx)  % the loop
      
x=xx{qq};

len_x = length(x);
	fltr=[1 1 1]/3;

    x1=x(1); x2=x(len_x); % remember the first and last elements. Irrelevant if not filtering.
    
  	for jj=1:filt  % run the averaging filter "filt" times (possibly zero)
	c=conv(fltr,x);
	x=c(2:len_x+1);
	x(1)=x1;  
        x(len_x)=x2; 
	end

lmval=[];
indd=[];
i=2;		% start at second data point in time series

    while i < len_x-1,
	if x(i) < x(i-1)
	   if x(i) < x(i+1)	% definite min
lmval =[lmval x(i)];
indd = [ indd i];

	   elseif x(i)==x(i+1)&x(i)==x(i+2)	% 'long' flat spot
%lmval =[lmval x(i)];	%1   comment these two lines for strict case 
%indd = [ indd i];	%2 when only  definite min included
i = i + 2;  		% skip 2 points

	   elseif x(i)==x(i+1)	% 'short' flat spot
%lmval =[lmval x(i)];	%1   comment these two lines for strict case
%indd = [ indd i];	%2 when only  definite min included
i = i + 1;		% skip one point
	   end
	end
	i = i + 1;
    end

if filt>0 & ~isempty(indd)
	if (indd(1)<= 3)|(indd(length(indd))+2>length(xx)) 
	   rng=1;	%check if index too close to the edge
	else rng=2;
	end

	   for ii=1:length(indd) 
		[val(ii) iind(ii)] = min(xx(indd(ii) -rng:indd(ii) +rng));
		iind(ii)=indd(ii) + iind(ii)  -rng-1;
	   end
  indd=iind; lmval=val;
else
end

lmvals{qq} = lmval;
indds{qq} = indd;
end

