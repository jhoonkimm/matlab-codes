function [SIcomp, SImag] = hodographsymmetryindex(v);

% HODOGRAPHSYMMEYRTINDEX calculate gait symmetry of hodographs
%  Based on method by Peter Adamczyk (2009) using
%  cross correlation.

% Peter Adamczyk

npstride = size(v,1);

%% Overall Symmetry using XCORR on all three components of FRCvelCOM
%% xcorr is actually called with the DEVIATIONS of FRCvelCOM from average walking speed 
[cx,lags] = xcorr(v(:,1)-mean(v(:,1)),floor(3/4*npstride),'unbiased');
[cy,lags] = xcorr(v(:,2)-mean(v(:,2)),floor(3/4*npstride),'unbiased');
[cz,lags] = xcorr(v(:,3)-mean(v(:,3)),floor(3/4*npstride),'unbiased');  % "floor(3/4*npstride)" gives the maximum lag. The "unbiased" option normalizes by the number of points overlapping for each lag, but it has artifacts at the ends where there aren't very many.  
c = [cy cz]; %2D  %%% c = [abs(cx) cy cz]; %3D  %%%  use either 2 or 3 components
%
%% Overall Symmetry Index using By-Component Normalization (mean of correlation coefficients)
ctot = mean(c./(ones(size(cz))*c(lags==0,:)),2);  % dividing each by the value at zero lag normalizes the xcorr to a max value of 1 (becomes an (unbiased) correlation coefficient)
[ctotmax,lagctotmax] = lmax(ctot);
ctotmax = ctotmax(abs(lags(lagctotmax))>10);  % look only at values with more than 10 sample lag 
lagctotmax = lagctotmax(abs(lags(lagctotmax))>10); % This and next line get rid of middle peak (if necessary)
[SIcomp,lagcompind] = max(ctotmax);
SIcomp;
lagcomp = abs(lags(lagctotmax(lagcompind)));
%
%% Overall Symmetry Index using Magnitude Normalization (normalizing factor is the Sum-Squared-Total of the Magnitude ||v-vavg|| )  
ctot = sum(c,2)/sum(c(lags==0,:));  % "sum(c(lags==0,:))" gives the SST of the magnitude of the velocity deviations from vavg
[ctotmax,lagctotmax] = lmax(ctot);
ctotmax = ctotmax(abs(lags(lagctotmax))>10);  % look only at values with more than 10 sample lag 
lagctotmax = lagctotmax(abs(lags(lagctotmax))>10); % This and next line get rid of middle peak (if necessary)
[SImag,lagmagind] = max(ctotmax);
SImag;
lagmag = abs(lags(lagctotmax(lagmagind)));
