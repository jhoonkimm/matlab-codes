function [p_A,p_B,p_AxB] = repAnova211(X,reps)
%REPANOVA211 2 factor anova with between-subject first factor and
%              within-subject (i.e., repeated) second factor
%
%   REPANOVA211(X,REPS) performs a balanced two-way ANOVA for comparing the
%   means of two or more columns and two or more rows of the sample in X.
%   The data in different columns represent changes in one factor. The data
%   in different rows represent changes in the other factor. If there is 
%   more than one observation per row-column pair, then the argument, REPS,
%   indicates the number of observations per "cell". A cell contains REPS 
%   number of rows.
%
%   First (independent) factor has levels down the rows, second (repeated) factor
%   has levels along the columns
%
%   For example, if REPS = 3, then each cell contains 3 rows and the total
%   number of rows must be a multiple of 3. If X has 12 rows, and REPS = 3,
%   then the "row" factor has 4 levels (3*4 = 12). The second level of the 
%   row factor goes from rows 4 to 6.  
%
% based on MATLAB's ANOVA2 function and 
% Bortz, J. Statistik fuer Naturwissenschaftler, 4th ed. Springer. pp. 307
%
% includes df correction for non-circular data
%
% Maximilian Riesenhuber 10/21/00

if nargin == 1,
  reps = 1;
end

if(floor(size(X,1)/reps)~=size(X,1)/reps)
  error('The number of rows must be a multiple of reps.');    
end

p=size(X,1)/reps;			% row factor A
q=size(X,2);		% column factor B
n=reps;					%  # of subjects

% put data into 3D matrix for convenience, last index is subject #
xijm=zeros(p,q,reps);
for m=1:reps
  xijm(:,:,m)=X(m+(0:reps:(p-1)*reps),:);
end

% calculate various values
G=sum(sum(sum(xijm)));				% total sum
Pim=squeeze(sum(xijm,2));
Ai=squeeze(sum(sum(xijm,2),3));
Bj=squeeze(sum(sum(xijm,1),3));
ABij=squeeze(sum(xijm,3));

% calculate the magic 8 terms
t1=G^2/(n*p*q);
t2=sum(sum(sum(xijm.^2)));
t3=sum(Ai.^2)/(q*n);
t4=sum(Bj.^2)/(p*n);
t5=sum(sum(ABij.^2))/n;
t6=sum(sum(Pim.^2))/q;

SS_A=t3-t1;
SS_inC=t6-t3;
SS_btwS=t6-t1;
SS_B=t4-t1;
SS_AxB=t5-t4-t3+t1;
SS_BxS=t2-t5-t6+t3;
SS_inS=t2-t6;
SS_tot=t2-t1;

df_A=p-1;
df_inC=p*(n-1);
df_btwS=p*n-1;
df_B=q-1;
df_AxB=(p-1)*(q-1);
df_BxS=p*(q-1)*(n-1);
df_inS=p*n*(q-1);
df_tot=p*q*n-1;

sig_A=SS_A/df_A;
sig_inC=SS_inC/df_inC;
sig_btwS=SS_btwS/df_btwS;
sig_B=SS_B/df_B;
sig_AxB=SS_AxB/df_AxB;
sig_BxS=SS_BxS/df_BxS;
sig_inS=SS_inS/df_inS;
sig_tot=SS_tot/df_tot;

% correct degrees of freedom, only for repeated factor
epsA=correctDOF(xijm,1);

df_A=floor(epsA*df_A);
df_inC=floor(epsA*df_inC);
df_AxB=floor(epsA*df_AxB);
df_BxS=floor(epsA*df_BxS);

if(sig_inS~=0)
  F_A=sig_A/sig_inC;
  p_A= 1- fcdf(F_A,df_A*epsA,df_inC*epsA);
elseif(sig_A==0)
  F_A=0;
  p_A=1;
else
  F_A = Inf;
  p_A=0;
end

if(sig_BxS~=0)
  F_B=sig_B/sig_BxS;
  p_B= 1- fcdf(F_B,df_B,df_BxS);
elseif(sig_B==0)
  F_B=0;
  p_B=1;
else
  F_B = Inf;
  p_B=0;
end

if(sig_AxB~=0)
  F_AxB=sig_AxB/sig_BxS;
  p_AxB= 1- fcdf(F_AxB,df_AxB*epsA,df_BxS*epsA);
elseif(sig_AxB==0)
  F_AxB=0;
  p_AxB=1;
else
  F_AxB = Inf;
  p_AxB=0;
end

if(isnan(p_A))				% this case happens when correction
  p_A=1;				% of df leads to 0 df
end
if(isnan(p_B))
  p_B=1;
end
if(isnan(p_AxB))
  p_AxB=1;
end

numrows = 8;
Table=zeros(numrows,4);   %Formatting for ANOVA Table printout with interactions.
Table(:,1)=[ SS_A SS_inC SS_btwS SS_B SS_AxB SS_BxS SS_inS SS_tot]';
Table(:,2)=[ df_A df_inC df_btwS df_B df_AxB df_BxS df_inS df_tot]';
Table(:,3)=[ sig_A sig_inC sig_btwS sig_B sig_AxB sig_BxS sig_inS sig_tot]';
Table(:,4)=[ F_A F_B F_AxB Inf Inf Inf Inf Inf]';
  
figh = figure('pos', [50 350 500 300]);
z0 = 0.00;
y0 = 0.85;
dz = 0.15;
dy =0.06;
text(z0+0.40,y0,'ANOVA Table')
axis off

z=z0+dz;
y=y0-3*dy/2;
if(epsA~=1)
  colheads = ['Source       ';'         SS  ';'   df (corr.)';'       MS    ';...
	'          F  '];
else
  colheads = ['Source       ';'         SS  ';'          df ';'       MS    ';...
'          F  '];
end

for i=1:5
  text(z,y,colheads(i,:))
  z = z + dz;
end
drawnow
rowheads = ['A          ';'inC        ';'btwS       ';'B          ';'AxB        ';'BxS        ';'inS        ';'tot        '];
for i=1:numrows
  y = y0-(i+1.5)*dy;
  z = z0 + dz;
  h = text(z,y,rowheads(i,:));
  
  z = z + dz;
  for j=1:4
    if (Table(i,j) ~= Inf) & j ~=2,
      h = text(z,y,sprintf('%11.4g',Table(i,j)));
      z = z + dz;
    elseif j==2,
      z = z + dz/4;
      h = text(z,y,sprintf('%7.5g',Table(i,j)));
      z = z + 3*dz/4;
    end
  end
  drawnow
end

function eps=correctDOF(x,factor)
  if(factor~=1 & factor ~=2)
    error('factor has to be 1 or 2!');
  end
  
% calculate covariance matrix
p=size(x,1);			% column factor A
q=size(x,2);		        % row factor B
n=size(x,3);				% number of subjects in each cell

% get covariance matrices for each level of the factor and average
for l=1:size(x,factor)
  if(factor==1)
    s(:,:,l)=cov(squeeze(x(l,:,:))');
  else
    s(:,:,l)=cov(squeeze(x(:,l,:))');
  end
end
sigIJ=squeeze(mean(s,3));
sigI=mean(sigIJ);
meanSigII=mean(diag(sigIJ));
meanSig=mean(mean(sigIJ));

% calculate epsilon
e=1/(p-1)*(p^2*(meanSigII-meanSig)^2)/(sum(sum(sigIJ.^2))-2*p*sum(sigI.^2)+p^2*meanSig^2);
eps=1;					% default: no correction

if(e<0.75)
  eps=e;
elseif(e<1.0)				% epsilon correction of Huynh and Feldt for 2 factors
  eTilde=(p*n*(q-1)*e-2)/((q-1)*(p*n-p-(q-1)*e));
  if(eTilde<1)
    eps=eTilde;
  end
end
