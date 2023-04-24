% Plots the limits p2L p2U p1L
clear all
close all
clc

L = 1;
mu1 = 7/4;
mu2 = mu1/5;
theta = 5;
g = 1/4;
c = 2/3;

pmax = (2*c / mu2) * ((theta - 1)/exp(1))^(theta - 1);
xLim = 2*pmax;
yLim = 2*pmax;

Fbar = @(p1) (exp(-0.05*p1));

p2L = @(p1, p2) p2 - (g*p1 .* max( (L/mu1)./Fbar(p1) - 1, 0));

p2U = @(p1, p2) (p2 + g*p1 - sqrt(g * p1 * pmax - g*(1-2*g)* p1.*p1));

p1L = @(p1, p2) (p1 + p2/g - 1/(g*sqrt(1-g)) * sqrt((1-2*g)*p2.*p2 + g*p2*pmax));


figure(1);
xlim([0 xLim])
ylim([0 yLim])
box on
hold on

% Shade the area (this is not a recommended method!)

for pp1 = linspace(0, pmax/(1-g), 60)
    for pp2 = linspace(0, pmax, 50)
        if (p1L(pp1, pp2) <0) && (p2U(pp1, pp2) <0) && (p2L(pp1, pp2) > 0)
            plot(pp1, pp2, 'k.')
        end
    end
end

hFig = gca;
hFig.ColorOrderIndex = 1;

% Draw the lines

hp1L = fimplicit(p1L, [0 xLim 0 yLim], 'Linewidth', 3);
hp2U = fimplicit(p2U, [0 xLim 0 yLim], 'Linewidth', 3);
hp2L = fimplicit(p2L, [0 xLim 0 yLim], 'Linewidth', 4);


fimplicit(@(p1, p2) p2 - pmax, [0 xLim 0 yLim], 'k-.', 'Linewidth', 1), 
fimplicit(@(p1, p2) p1 - (1/(1-g)) * pmax, [0 xLim 0 yLim], 'k-.', 'Linewidth', 1), 


% Label the plots.

xlabel('$p_1$', 'Interpreter', 'Latex', 'Fontsize', 24),
ylabel('$p_2$', 'Interpreter', 'Latex', 'Fontsize', 24),
lgnd = legend([hp1L, hp2U, hp2L], {'$p_1^{\sf L}$','$p_2^{\sf U}$','$p_2^{\sf L}$'}, ...
    'Interpreter', 'Latex', ...
    'Fontsize', 24, ...
    'Orientation','horizontal', ...
    'Position',[0.142806535139881 0.765151515151515 0.330239286962545 0.141098459636744], ...
    'Color', [0.9 0.9 0.9], ...
    'EdgeColor',[1 1 1]);

text(pmax*2/3, pmax*10/9, '$p_{\max}$', ...
    'Interpreter', 'Latex', 'Fontsize', 24)

text(1.05*pmax/(1-g), pmax*3/2, '$\frac{1}{1-\gamma}p_{\max}$', ...
    'Interpreter', 'Latex', 'Fontsize', 24)

hold off

saveas(gcf,'../Figures/plotSearch1.png')
% system('epstopdf ../Figures/plotSearch2.eps')