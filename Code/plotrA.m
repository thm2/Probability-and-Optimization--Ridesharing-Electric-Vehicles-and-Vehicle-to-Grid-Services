% Plots r_A

L = 1;
mu1 = 1.2;
mu2 = 0.1*mu1;
theta = 5;
g = 0.25;
c = 10;

pUpper = 2*c;

rA = @(p1, p2) ((L * g*p1 ./ (g*p1 + p2)) .* (1-g) .* p1 ...
    - (L * p2 ./ (g*p1 + p2)) .* p2 ...
    + c * (1 - 2*gammainc(( (L/mu2) * p2 ./ (g*p1 + p2)) , theta, 'upper'))) ;

meshDensity = 35;
    % For the final version, increase it to 100.

h = fsurf(rA, [0 pUpper 0 pUpper],'ShowContours','on','MeshDensity',meshDensity)
view(-19, 31),
colormap default

% Uncomment the next four lines for the final version of the figure.

camlight(110,70)
brighten(0.01)
h.EdgeColor = 'none';
h.AmbientStrength = 0.6;


a = gca;
a.TickLabelInterpreter = 'latex';
a.Box = 'on';
a.BoxStyle = 'full';

xlabel('$p_1$', 'Interpreter', 'Latex', 'Fontsize', 24)
ylabel('$p_2$', 'Interpreter', 'Latex', 'Fontsize', 24),
zlabel('$r_{\mathbf A}(p_1, p_2)$', 'Interpreter', 'Latex', 'Fontsize', 24)

saveas(gcf,'../Figures/plotrA.png')
% system('epstopdf ../Figures/plotrA.eps')
