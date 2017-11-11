function [patchhandles, linehandles] = errpatch2dcovellipse(x,y,arrdims,mndim,varargin)

% ERRPATCH2COVELLIPSE plot error bounds around a curve in 2D
%
% Function plots the error bounds around a curve in 2 dimensions 
% by drawing the Covariance Ellipse at each point:  
%
% Function call: errpatch2dcovellipse(X, Y, arrdims, mndim, [Color], [LineSpec], ['closed'])
%
%   X, Y: arrays holding the point-by-point X and Y data for several samples
%     of a curve, which will be averaged along "mndim". X and Y can have no
%     more than 3 non-singleton dimensions
%   arrdims: dimensions of the arrays to input to a simple 2-D "plot" command
%     (2 of them) [1 2]
%   mndim: dimension of the arrays X and Y should we compute mean
%     and covariance. [3]
%
% VARARGIN:
%   Color: rgb triples in 2nd or 3rd dim, stacked in 2nd dim if more than one
%   LineSpec: a line specification per the PLOT command
%   'Closed': input this string to have the PATCH that is produced close
%      perfectly on the beginning and end of the curve. 

% Peter Adamczyk

clr = [.75 .75 .75];
lspc = '-';
closed = 0;

for ii = 1:nargin-4
    if ~isempty(varargin{ii})
        switch ii
            case 1
                clr = varargin{ii};
            case 2
                lspc = varargin{ii};
            case 3
                closed = 1;
        end
    end
end

% Need to know number of dims in order to reshape matrices properly. 
nd = ndims(x);

if nd ~= 3 
    errmsg('Array must have 3 dimensions (PGA)')
end

nshift = mndim;  % shift dimensions to make sure the "across this" dimension is LAST
shiftx = permute(x,[arrdims mndim]); 
shifty = permute(y,[arrdims mndim]); 

% Compute mean of X and Y, 
mnx = mean(shiftx,nd); 
mny = mean(shifty,nd);  
% Compute covariance matrix between X and Y. 
for ii = 1:size(shiftx,1)
    for jj = 1:size(shiftx,2)
        covariancematrices{ii,jj} = cov(squeeze(shiftx(ii,jj,:)),squeeze(shifty(ii,jj,:)));
    end
end

% Find the Normal to the Curve at each point
% First, the Average Derivative at each point
xdiff = diff(mnx); xdiff = mean(cat(3,[xdiff(1,:);xdiff],[xdiff;xdiff(end,:)]),3); % mean difference at each point, including Endpoints. 
ydiff = diff(mny); ydiff = mean(cat(3,[ydiff(1,:);ydiff],[ydiff;ydiff(end,:)]),3);
% now the Normal Vector (-dX and dY form vectors along matrix dimension 3)
normalvec = cat(3,-ydiff,xdiff); % "Inward" Normal Vector
% % % Normalize the Normal
% % unitnormalvec = normalvec./repmat(vecnorm(normalvec,3),[1 1 2]);

% % % Break the Normal Vector into cells
% % unitnormalveccell = num2cell(unitnormalvec,3);
normalveccell = num2cell(normalvec,3);


%%%%% PETER's Method Obsolete as of 2008-09-23. The Method was supposed to 
%%%%% pick out the point on each ellipse that had a tangent parallel to the
%%%%% curve and use only that. The algorithm did this, but the
%%%%% covariance ellipses spilled over the chosen points at sharp turns,
%%%%% probably because of the difference between actual and estimated
%%%%% Tangents and Normals.
% % % 
% % % use Peter's method to produce a vector to the point on the Covariance
% % % ellipse that is farthest from the curve point in the normal direction   
% % 
% % % %% Below is the primitive for the algorithm:
% % % %% See Peter's Notes of 9/19/2008-9/21/2008
% % % [v,d] = eig(covmat);
% % % vsqrtd = v*sqrt(d)
% % % invvsqrtd = inv(vsqrtd)
% % % flop = [0 -1;1 0]
% % % n = [1 1];  % EXAMPLE n vector. It 
% % % K = vsqrtd*flop*invvsqrtd 
% % % nk = n*K
% % % rdir = cross([nk 0], [0 0 1])  % transformed "Inward" Normal >cross< +Z "binormal" vector gives Tangent to ellipse  
% % % R = vecnorm(invvsqrtd*rdir(1:2)')  
% % % rellipse = rdir(1:2)*1/R
% % % 
% % % outerpoint = [0 0]+rellipse
% % % innerpoint = [0 0]+ -1*rellipse
% % % blah = [cos(t); sin(t)];
% % % ellipsepts = vsqrtd*blah;
% % % figure; plot(ellipsepts(1,:),ellipsepts(2,:)); hold on; axis equal
% % % plot([outerpoint(1);0;innerpoint(1)], [outerpoint(2); 0; innerpoint(2)],'ko')
% % % outer2 = outerpoint'+K*rellipse';
% % % inner2 = innerpoint'+-1*K*rellipse';
% % % plot([outerpoint(1);outer2(1)],[outerpoint(2);outer2(2)],'gx-')
% % % plot([innerpoint(1);inner2(1)],[innerpoint(2);inner2(2)],'rx-')
% % 
% % % Algorithm Implemented using Cell Arrays
% % [v,d] = cellfun(@eig,covariancematrices,'UniformOutput',0);
% % vsqrtd = cellfun(@(vv,dd) vv*sqrt(dd),v,d,'UniformOutput',0);
% % invvsqrtd = cellfun(@inv,vsqrtd,'UniformOutput',0);
% % flop = [0 -1 ; 1 0];
% % K = cellfun(@(a,ainv) a*flop*ainv, vsqrtd, invvsqrtd, 'UniformOutput',0);
% % normalveccellK = cellfun(@(n,k) reshape(n, [1 2])*k, normalveccell, K,'UniformOutput',0);
% % rdir = cellfun(@(nk) cross([reshape(nk,[1 2]) 0],[0 0 1]), normalveccellK, 'UniformOutput',0);
% % rdir = cellfun(@(blah) blah(1:2), rdir, 'UniformOutput',0);
% % R = cellfun(@(iinvvsqrtd, rrdir) vecnorm(iinvvsqrtd*rrdir'), invvsqrtd, rdir, 'UniformOutput',0); 
% % rellipse = cellfun(@(rrdir, RR) rrdir*1/RR, rdir, R, 'UniformOutput',0);
% % 
% % % Create vectors for Normal and NegativeNormal offsets from the original curve
% % xinwardnormal = cellfun(@(x,r) x + r(1), num2cell(mnx), rellipse, 'UniformOutput',1);
% % xoutwardnormal = cellfun(@(x,r) x + -1*r(1), num2cell(mnx), rellipse, 'UniformOutput',1);
% % yinwardnormal = cellfun(@(x,r) x + r(2), num2cell(mny), rellipse, 'UniformOutput',1);
% % youtwardnormal = cellfun(@(x,r) x + -1*r(2), num2cell(mny), rellipse, 'UniformOutput',1);
% % 
% % % plot([xinwardnormal; xoutwardnormal],[yinwardnormal;youtwardnormal],'k.')
% % 
% % if closed
% %     xfill = [xoutwardnormal([1:end,1],:); xinwardnormal([1,end:-1:1],:)];
% %     yfill = [youtwardnormal([1:end,1],:); yinwardnormal([1,end:-1:1],:)];
% % else
% %     xfill = [xoutwardnormal;xinwardnormal(end:-1:1,:)];
% %     yfill = [youtwardnormal;yinwardnormal(end:-1:1,:)];
% % end
% % 
% % out = fill(xfill, yfill, reshape(clr,1,[],3));
% % set(out,'EdgeColor','none')


%% New Method using Part of Peter's method, but only to draw the Covariance 
%% Ellipse directly. 
% Algorithm Implemented using Cell Arrays
[v,d] = cellfun(@eig,covariancematrices,'UniformOutput',0);
vsqrtd = cellfun(@(vv,dd) vv*sqrt(dd),v,d,'UniformOutput',0);
t = [-pi:pi/99:pi];
rcirc = [cos(t); sin(t)];
% % Plot the ellipses
hold on;
ellipsehandles = cellfun(@(M,mx,my) fill(M(1,:)*rcirc + mx, M(2,:)*rcirc + my, reshape(clr,1,[],3)), vsqrtd, num2cell(mnx), num2cell(mny), 'UniformOutput', 1);
set(ellipsehandles,'EdgeColor','none')

ncurves = size(ellipsehandles, 2);
clr = linspace(1,0.5,ncurves)'*clr;
for ii = 1:ncurves
    set(ellipsehandles(:,ii),'FaceColor',clr(ii,:));

    x0 = get(ellipsehandles(:,ii),{'XData','YData'}); % get a cell array of data, each cell a vector for one ellipse
    [x,y] = deal(x0(end,1),x0(end,2));
    for jj = 1:length(ellipsehandles)
        [x,y] = polybool('union',x,y,x0(jj,1),x0(jj,2)); % compute the union of all the ellipses
    end
    [f,v] = poly2fv(x,y);
    errpatchhandles(ii) = patch('Faces',f,'Vertices',v,'FaceColor',clr(ii,:),'EdgeColor','none');
    delete(ellipsehandles(:,ii));
end


out2 = plot(mnx,mny,lspc);
% patchhandles = ellipsehandles; 
patchhandles = errpatchhandles; 
linehandles = out2;
hold off
% keyboard

end
