function plotEllipse(stepWidth, stepLength)
% two input vectors are the x and y locations of the data

% some other stuffs
stanDev = 1;                 % specifies # of stanDevs ellipse will cover
conf = 2*normcdf(stanDev)-1;     % covers around 95% of population
scale = chi2inv(conf,2);     % inverse chi-squared with dof=#dimensions

Cov = cov(stepWidth, stepLength)*scale;
[V D] = eig(Cov);

t = linspace(0,2*pi,100);
e = [cos(t) ; sin(t)];        % unit circle
VV = V*sqrt(D);               % scale eigenvectors
e = bsxfun(@plus, VV*e, 0);   % project circle back to orig space

% plot cov and major/minor axes
plot(e(1,:), e(2,:), 'Color','k','LineWidth',2);
xlabel('Width'); ylabel('Length')

% Specify axes
xmin = -100; xmax = 100;
ymin = -100; ymax = 100;

axis equal
axis([xmin xmax ymin ymax])

end