clear all;
close all;

% simulation parameters
g1=0.25; m01=3; mu2=3; theta=10; c=15; lambda=1;
pmax=2*c/(mu2*gamma(theta))*exp(1-theta)*((theta-1)^(theta-1));

% parameters for gamma distribution used in stability constraint
ag=2; bg=1;

% axis for p_1<1.33 pmax; % bottom branch of proposition 2
p1_start=0.0; dt1=0.005; p1_end=1/(1-g1)*pmax-dt1; p1=linspace(p1_start,p1_end,(p1_end-p1_start)/dt1);

% CURVES %
% p_2^L
curve1=max(g1.*p1.*(lambda./(m01*(1-gamcdf(p1,ag,bg)))-1),0);

% p_2^U
curve2=-g1.*p1+sqrt(g1.*p1*pmax-g1*(1-2*g1).*p1.^2);

% PLOTS for bottom branch %
plot(p1,curve2,'LineWidth',1.5)
hold on
plot(p1,curve1,'LineWidth',1.5)
hold on

% this paints everything inbetween curve2 and curve1
x2 = [p1, fliplr(p1)];
inBetween = [curve1, fliplr(curve2)];
fill(x2, inBetween, 'g');

% this only paints the difference curve2-curve1
% idx2 = curve2 > curve1 ;
% xb = p1(idx2) ;
% y1b = curve1(idx2) ;
% y2b = curve2(idx2) ;
% patch([xb fliplr(xb)],[y1b fliplr(y2b)], 'g')


% axis for p_1>1.33 pmax; top branch of proposition 2
p_start=1/(1-g1)*pmax;p_end=2*p_start; p=p_start:dt1:p_end;

% PLOT for top branch %
plot(p,max(g1.*p.*(lambda./(m01*(1-gamcdf(p,ag,bg)))-1),0),'g','LineWidth',1.5);

a = gca;
a.TickLabelInterpreter = 'latex';
% a.Box = 'on';
% a.BoxStyle = 'full';

xlabel('$p_1$', 'Interpreter', 'Latex', 'Fontsize', 24)
ylabel('$p_1$', 'Interpreter', 'Latex', 'Fontsize', 24)
title('Search space', 'Interpreter', 'Latex', 'Fontsize', 24)