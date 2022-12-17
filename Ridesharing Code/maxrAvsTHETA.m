clear all
close all

RAMAXFINAL=[];
RAMAX1=[];

%variables
g1=0.25;
lambda=20;
mu1=10; %mu1 of tex file
mu2=2;
c=15;
theta=5;

theta=[1;2;3;5;10;15;20];


maxcfp2=(g1-sqrt(-g1*(-1+g1)/(-1+2*g1)^2)+2*g1*sqrt(-g1*(-1+g1)/(-1+2*g1)^2))/(2*(-1+2*g1));
maxcfp1=1/(1-g1);

% gamma parameters
ag=2; bg=1; E=ag*bg;

% axes - first max
% p1_start=0.01; dt1=0.05; p1_end=maxcfp1*pmax-dt1; p1=linspace(p1_start,p1_end,(p1_end-p1_start)/dt1);
% p2_start=0.01; dt2=0.05; p2_end=maxcfp2*pmax; p2=linspace(p2_start,p2_end,(p2_end-p2_start)/dt2);
% l_start=1; dt3=20; l_end=700; lambda=linspace(l_start,l_end,(l_end-l_start)/dt3);

% arrival rates
lambda1=@(p1,p2,lambda)lambda*(g1*p1/(g1*p1+p2));
lambda2=@(p1,p2,lambda)lambda*(p2/(g1*p1+p2));

lamm=0; mn=1;

clc
RMAX=[]; P22=[]; P11=[];
for lambda=8:2:12
    RAMAX1=[]; ERR=[]; 
lamm=0;  
for theta_incr=1:length(theta) %
    pmax=2*c/(mu2*gamma(theta(theta_incr)))*(theta(theta_incr)-1)^(theta(theta_incr)-1)*exp(1-theta(theta_incr));
    p1_start=0.01; dt1=0.02; p1_end=maxcfp1*pmax-dt1; p1=linspace(p1_start,p1_end,(p1_end-p1_start)/dt1);
p2_start=0.01; dt2=0.02; p2_end=maxcfp2*pmax; p2=linspace(p2_start,p2_end,(p2_end-p2_start)/dt2);
rA1=zeros(length(p1),length(p2),length(theta));
        lamm=lamm+1;

    for jj=1:length(p2)
        for ii=1:length(p1) %
           ER = c*(1-2*gammainc(lambda2(p1(ii),p2(jj),lambda)/mu2,theta(theta_incr),'upper'));
           rA1(ii,jj,lamm)=((1-g1)*lambda1(p1(ii),p2(jj),lambda)*p1(ii)+ ER -lambda2(p1(ii),p2(jj),lambda)*p2(jj))/mn;      
        end
    end

    RAMAX11=max(max(rA1(:,:,lamm)));
                display('k')

    RAMAX1=[RAMAX1;max(max(RAMAX11))];
end




% LLL=[]; LINEARINDEX=[]; P1=[]; P2=[]; RAMAX1=[];
% % here it extracts the optimal prices from TotExpProfitRate
% for LL=1:size(rA1,3);
%      [maxValue, linearIndexesOfMaxes] = max(max(rA1(:,:,LL)));
%     if min(min(isnan((rA1(:,:,LL)))))~=1 % only the non nan exp. rates (they have at least one max) are kept
%         [rowsOfMaxes, colsOfMaxes] = find(rA1(:,:,LL) == maxValue)
%         LINEARINDEX=[LINEARINDEX; sub2ind(size(rA1), rowsOfMaxes, colsOfMaxes, LL)];
%         P1=[P1; p1(rowsOfMaxes)];
%         P2=[P2; p2(colsOfMaxes)];
%         LLL=[LLL; LL];
%     end
%     % first potential maximum
%     rAmax1=maxValue;
%     RAMAX1=[RAMAX1; rAmax1];
% end




%%%% Seconnd search space

PMAX2=[];
RAMAX2=[];

for theta_incr=1:length(theta) %

eta=0;
etaprime=mu1/g1*((1-g1)*mu1*E+c-eta);
psi=min(lambda/mu1,1);

if etaprime<0
    p1bar=gaminv(1-psi,ag,bg);
else
    p1bar=gaminv(1-psi/2,ag,bg)*max(1,1/(4*etaprime)*mu1^2*psi^2);
end


% gradient ascend

% gamma cdf and derivative
F=@(p1)gamcdf(p1,ag,bg); dF=@(p1)p1^(ag-1)*exp(-p1/bg)/(gamma(ag)*bg^ag);

% gradient (derivative) of r_A(p1,p2L(p1))
f=@(p1,lambda,theta,c)(g1+(g1*lambda)/(mu1*(-1+F(p1)))+(-1+g1)*mu1*(-1+F(p1))...
    -g1*F(p1)-2*g1*p1*dF(p1)+(-1+g1)*mu1*p1*dF(p1)...
    +(g1*p1*(-lambda+mu1+mu1*(-2+F(p1))*F(p1))*dF(p1))/(mu1*(-1+F(p1))^2))...
    +(2*c*exp(-(lambda-mu1+mu1*F(p1))/(mu2))*mu1*((lambda-mu1+mu1*F(p1))/(mu2))^(theta-1)*dF(p1))/(mu2*gamma(theta));

x_0 = 0.0; %Initial guess
alpha = @(k)1/k; %Step size
iteration_max = 1000;
tolerance = 10e-8;

%Initiliazation
iter = 1; x_new=0;

while true 
    grad = f(x_0(1),lambda,theta(theta_incr),c);
    x_new = x_0 + alpha(iter) * grad; %New solution
    iter = iter + 1;
    
    if (abs(x_new - x_0) < tolerance)
        break
    end
    
    x_0 = x_new; %Update old solution
    
    if iter > iteration_max || x_0> p1bar
        break
    end
end

p1max2=x_0; p2max2=g1*x_0*max((lambda)/(mu1*(1-F(x_0)))-1,0);
rA2 = @(p1, p2, L, theta)((L * g1*p1./(g1*p1+p2)).*(1-g1).* p1 ...
    - (L * p2 ./ (g1*p1 + p2)) .* p2 ...
    + c * (1 - 2*gammainc(( (L/mu2) * p2 ./ (g1*p1 + p2)) , theta, 'upper')));
rAmax2=rA2(p1max2,p2max2,lambda,theta(theta_incr));

PMAX2=[PMAX2; p2max2 p1max2];
RAMAX2=[RAMAX2;rAmax2];
end
RMAX=[RMAX;max(RAMAX1,RAMAX2)];
% P11=[P11; P1]
% P22=[P22;P2]
end


figure
% plot(c,RMAX(1:length(c)),c,RMAX(length(c)+1:2*length(c)),c,RMAX(2*length(c)+1:3*length(c)),'LineWidth',3,strcat('$\lambda =$',num2str(8)))
plot(theta,RMAX(1:length(theta)),'r','LineWidth',3,'DisplayName',strcat('$\lambda =$',num2str(8)))
hold on
plot(theta,RMAX(length(theta)+1:2*length(theta)),'b','LineWidth',3, 'DisplayName', strcat('$\lambda =$',num2str(10)))
hold on
plot(theta,RMAX(2*length(theta)+1:3*length(theta)),'k','LineWidth',3,'DisplayName',strcat('$\lambda =$',num2str(12)))
hold on
% plot(c,RMAX(3*length(c)+1:4*length(c)),'LineWidth',3)
xlabel('$\theta$', 'Interpreter', 'Latex', 'Fontsize', 20)
ylabel('$r_{\mathbf A}^*$', 'Interpreter', 'Latex', 'Fontsize', 20),
title('Optimal $r_{\bf A}$ versus threshold $\theta$', 'Interpreter', 'Latex', 'Fontsize', 17)
% leg = legend(gca, 'show'), set(leg, 'fontsize',18, 'Interpreter','latex') % for dispname
grid on

saveas(gcf,'../Figures/maxrA_vs_theta.png')

% plot(theta,max(RAMAX1,RAMAX2),'LineWidth',3)
% xlabel('$\theta$', 'Interpreter', 'Latex', 'Fontsize', 24)
% ylabel('$r_{\mathbf A}(p_1, p_2)$', 'Interpreter', 'Latex', 'Fontsize', 24),
% title('Total expected reward rate versus contract threshold', 'Interpreter', 'Latex', 'Fontsize', 16)
% grid on
% saveas(gcf,'../Figures/maxrAvstheta.pdf','epsc')