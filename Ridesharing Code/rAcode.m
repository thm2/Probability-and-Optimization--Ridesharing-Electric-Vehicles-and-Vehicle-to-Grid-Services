%% FIRST CELL: VARYING LAMBDA
clear all
close all
clc


% axis parameters
p1_start=0.1; dt1=0.25; p1_end=6; p1=linspace(p1_start,p1_end,(p1_end-p1_start)/dt1);
p2_start=0.1; dt2=0.25; p2_end=6; p2=linspace(p2_start,p2_end,(p2_end-p2_start)/dt2);
l_start=1; dt3=20; l_end=700; lambda=linspace(l_start,l_end,(l_end-l_start)/dt3);
% lambda=[1;20;50;200;400;1000;1400;2000];

% simulation parameters
g1=0.25; m01=3000; mu2=3; theta=20; cc= 10000;

% arrival rates
lambda1=@(p1,p2,lambda)lambda*(g1*p1/(g1*p1+p2));
lambda2=@(p1,p2,lambda)lambda*(p2/(g1*p1+p2));
ag=2; bg=1; mu1=@(p1,p2)m01*(1-gamcdf(p1,ag,bg));

lamm=0; mn=1; TotExpProfitRate=zeros(length(p1),length(p2),length(lambda));


for lambda_incr=1:length(lambda) %
    lamm=lamm+1;
    
    for jj=1:length(p2)
        for ii=1:length(p1) %
            
            if lambda1(p1(ii),p2(jj),lambda(lambda_incr))<mu1(p1(ii),p2(jj)) % check transp. queue stability condition
                ER = cc*(1-2*gammainc(lambda2(p1(ii),p2(jj),lambda(lambda_incr))/mu2,theta,'upper'));
                TotExpProfitRate(ii,jj,lamm)=((1-g1)*lambda1(p1(ii),p2(jj),lambda(lambda_incr))*p1(ii)+ ER -lambda2(p1(ii),p2(jj),lambda(lambda_incr))*p2(jj))/mn;
                
            else % if transp. queue stability condition is not met, announce nothing
                TotExpProfitRate(ii,jj,lamm)=NaN;
                
            end
        end
    end
end

figure
for LL=1:size(TotExpProfitRate,3);
    surf(p1,p2,TotExpProfitRate(:,:,LL)','edgecolor', 'none');
    hold on
end
xlabel('p_1')
ylabel('p_2')
zlabel('revenue rate')
title('Total expected profit rate vs \lambda')


% code to plot the optimal prices versus arrival rate lambda
LLL=[]; LINEARINDEX=[]; P1=[]; P2=[];
% here it extracts the optimal prices from TotExpProfitRate
for LL=1:size(TotExpProfitRate,3);
    [maxValue, linearIndexesOfMaxes] = max(max(TotExpProfitRate(:,:,LL)));
    if min(min(isnan((TotExpProfitRate(:,:,LL)))))~=1 % only the non nan exp. rates (they have at least one max) are kept
        [rowsOfMaxes, colsOfMaxes] = find(TotExpProfitRate(:,:,LL) == maxValue);
        LINEARINDEX=[LINEARINDEX; sub2ind(size(TotExpProfitRate), rowsOfMaxes, colsOfMaxes, LL)];
        P1=[P1; p1(rowsOfMaxes)];
        P2=[P2; p2(colsOfMaxes)];
        LLL=[LLL; LL];
    end
end

%%% plot - maximizing prices vs lambda
figure
scatter(lambda(LLL)/mn,P1,[],TotExpProfitRate(LINEARINDEX),'filled','LineWidth',1.8)
hold on
scatter(lambda(LLL)/mn,P2,[],TotExpProfitRate(LINEARINDEX),'LineWidth',1.8,'MarkerEdgeColor',[0 .5 .5]) % ,
xlabel('\lambda')
ylabel('prices')
legend('p1*','p2*')
title('Prices for max profit vs \lambda')
grid on

%%% plot - Total profit rate at max prices vs lambda
figure
plot(lambda(LLL),TotExpProfitRate(LINEARINDEX),'-o','LineWidth',2)
xlabel('\lambda')
ylabel('revenue rate')
title('Total expected profit rate @ max prices vs \lambda')


%% SECOND CELL: VARYING THETA

clear all
close all
clc


% axis parameters
p1_start=0.1; dt1=0.25; p1_end=6; p1=linspace(p1_start,p1_end,(p1_end-p1_start)/dt1);
p2_start=0.1; dt2=0.25; p2_end=6; p2=linspace(p2_start,p2_end,(p2_end-p2_start)/dt2);
t_start=1; dt3=25; t_end=500; theta=linspace(t_start,t_end,(t_end-t_start)/dt3);
% theta=[1;20;60;100;200;500;750;1000];

% simulation parameters
g1=0.25; m01=3000; mu2=3; lambda=1000; cc= 10000;

% arrival rates
lambda1=@(p1,p2,lambda)lambda*(g1*p1/(g1*p1+p2));
lambda2=@(p1,p2,lambda)lambda*(p2/(g1*p1+p2));
ag=5; bg=1; mu1=@(p1,p2)m01*(1-gamcdf(p1,ag,bg));

lamm=0; mn=1; TotExpProfitRate=zeros(length(p1),length(p2),length(theta));


for theta_incr=1:length(theta) %
    lamm=lamm+1;
    
    for jj=1:length(p2)
        for ii=1:length(p1) %
            
            if lambda1(p1(ii),p2(jj),lambda)<mu1(p1(ii),p2(jj)) % check transp. queue stability condition
                ER = cc*(1-2*gammainc(lambda2(p1(ii),p2(jj),lambda)/mu2,theta(theta_incr),'upper'));
%                 ER = cc*(1-2*(igamma(theta(theta_incr),lambda*(p2(jj)/(g1*p1(ii)+p2(jj)))))/gamma(theta(theta_incr)));
                TotExpProfitRate(ii,jj,lamm)=((1-g1)*lambda1(p1(ii),p2(jj),lambda)*p1(ii)+ ER -lambda2(p1(ii),p2(jj),lambda)*p2(jj))/mn;
                
            else % if transp. queue stability condition is not met, announce nothing
                TotExpProfitRate(ii,jj,lamm)=NaN;
                
            end
        end
    end
end

figure
for LL=1:size(TotExpProfitRate,3);
    surf(p1,p2,TotExpProfitRate(:,:,6)');
    hold on
end
xlabel('p_1')
ylabel('p_2')
zlabel('revenue rate')
title('Total expected profit rate vs \theta = a representative theta')


% here it extracts the optimal prices from TotExpProfitRate
LLL=[]; LINEARINDEX=[]; P1=[]; P2=[];
for LL=1:size(TotExpProfitRate,3);
    [maxValue, linearIndexesOfMaxes] = max(max(TotExpProfitRate(:,:,LL)));
    if min(min(isnan((TotExpProfitRate(:,:,LL)))))~=1 % only the non nan exp. rates (they have at least one max) are kept
        [rowsOfMaxes, colsOfMaxes] = find(TotExpProfitRate(:,:,LL) == maxValue);
        LINEARINDEX=[LINEARINDEX; sub2ind(size(TotExpProfitRate), rowsOfMaxes, colsOfMaxes, LL)];
        P1=[P1; p1(rowsOfMaxes)];
        P2=[P2; p2(colsOfMaxes)];
        LLL=[LLL; LL];
    end
end

%%% plot - maximizing prices vs theta
figure
scatter(theta(LLL)/mn,P1,[],TotExpProfitRate(LINEARINDEX),'filled','LineWidth',1.8)
hold on
scatter(theta(LLL)/mn,P2,[],TotExpProfitRate(LINEARINDEX),'o','LineWidth',1.8,'MarkerEdgeColor',[0 .5 .5]) % ,
xlabel('\theta')
ylabel('prices')
legend('p1*','p2*')
title('Prices for max profit vs \theta')
grid on

%%% plot - Total profit rate at max prices vs theta
figure
plot(theta(LLL),TotExpProfitRate(LINEARINDEX),'-o','LineWidth',2)
xlabel('\theta')
ylabel('revenue rate')
title('Total expected profit rate @ max prices vs \theta')


%% THIRD CELL: VARYING MU01  (in mu1=mu01(1-Fv(p1);only affects constraint)


clear all
close all
clc


% axis parameters
p1_start=0.1; dt1=0.25; p1_end=6; p1=linspace(p1_start,p1_end,(p1_end-p1_start)/dt1);
p2_start=0.1; dt2=0.25; p2_end=6; p2=linspace(p2_start,p2_end,(p2_end-p2_start)/dt2);
m_start=1; dt3=30; m_end=2000; m01=linspace(m_start,m_end,(m_end-m_start)/dt3);
% mu01=[1;20;60;100;200;500;750;1000];

% simulation parameters
g1=0.25; mu2=3; lambda=1000; theta=20; cc=10000;

% arrival rates
lambda1=@(p1,p2,lambda)lambda*(g1*p1/(g1*p1+p2));
lambda2=@(p1,p2,lambda)lambda*(p2/(g1*p1+p2));
ag=5; bg=1; mu1=@(p1,p2)m01*(1-gamcdf(p1,ag,bg));

lamm=0; mn=1; TotExpProfitRate=zeros(length(p1),length(p2),length(m01));


for m_incr=1:length(m01) %
    lamm=lamm+1;
    
    for jj=1:length(p2)
        for ii=1:length(p1) %
            
            if lambda1(p1(ii),p2(jj),lambda) < m01(m_incr)*(1-gamcdf(p1(ii),ag,bg)) % check transp. queue stability condition
                ER = cc*(1-2*gammainc(lambda2(p1(ii),p2(jj),lambda)/mu2,theta,'upper'));
%                 ER = cc*(1-2*(igamma(theta(theta_incr),lambda*(p2(jj)/(g1*p1(ii)+p2(jj)))))/gamma(theta(theta_incr)));
                TotExpProfitRate(ii,jj,lamm)=((1-g1)*lambda1(p1(ii),p2(jj),lambda)*p1(ii)+ ER -lambda2(p1(ii),p2(jj),lambda)*p2(jj))/mn;
                
            else % if transp. queue stability condition is not met, announce nothing
                TotExpProfitRate(ii,jj,lamm)=NaN;
                
            end
        end
    end
end

figure
for LL=1:size(TotExpProfitRate,3);
    surf(p1,p2,TotExpProfitRate(:,:,LL)');
    hold on
end
xlabel('p_1')
ylabel('p_2')
zlabel('revenue rate')
title('Total expected profit rate vs \mu_{01}')


% here it extracts the optimal prices from TotExpProfitRate
LLL=[]; LINEARINDEX=[]; P1=[]; P2=[];
for LL=1:size(TotExpProfitRate,3);
    [maxValue, linearIndexesOfMaxes] = max(max(TotExpProfitRate(:,:,LL)));
    if min(min(isnan((TotExpProfitRate(:,:,LL)))))~=1 % only the non nan exp. rates (they have at least one max) are kept
        [rowsOfMaxes, colsOfMaxes] = find(TotExpProfitRate(:,:,LL) == maxValue);
        LINEARINDEX=[LINEARINDEX; sub2ind(size(TotExpProfitRate), rowsOfMaxes, colsOfMaxes, LL)];
        P1=[P1; p1(rowsOfMaxes)];
        P2=[P2; p2(colsOfMaxes)];
        LLL=[LLL; LL];
    end
end

%%% plot - maximizing prices vs mu01
figure
scatter(m01(LLL)/mn,P1,[],TotExpProfitRate(LINEARINDEX),'filled','LineWidth',1.8)
hold on
scatter(m01(LLL)/mn,P2,[],TotExpProfitRate(LINEARINDEX),'o','LineWidth',1.8,'MarkerEdgeColor',[0 .5 .5]) % ,
xlabel('\mu_{01}')
ylabel('prices')
legend('p1*','p2*')
title('Prices for max profit vs \mu_{01}')
grid on

%%% plot - Total profit rate at max prices vs lambda
figure
plot(m01(LLL),TotExpProfitRate(LINEARINDEX),'-o','LineWidth',2)
xlabel('\mu_{01}')
ylabel('revenue rate')
title('Total expected profit rate @ max prices vs \mu_{01}')


%% FOURTH CELL: VARYING cc


clear all
close all
clc


% axis parameters
p1_start=0.1; dt1=0.25; p1_end=6; p1=linspace(p1_start,p1_end,(p1_end-p1_start)/dt1);
p2_start=0.1; dt2=0.25; p2_end=6; p2=linspace(p2_start,p2_end,(p2_end-p2_start)/dt2);
c_start=1; dt3=30; c_end=2000; cc=linspace(c_start,c_end,(c_end-c_start)/dt3);
% cc=[1;20;60;100;200;500;750;1000];

% simulation parameters
g1=0.25; mu2=3; lambda=1000; theta=20; m01=3000;

% arrival rates
lambda1=@(p1,p2,lambda)lambda*(g1*p1/(g1*p1+p2));
lambda2=@(p1,p2,lambda)lambda*(p2/(g1*p1+p2));
ag=5; bg=1; mu1=@(p1,p2)m01*(1-gamcdf(p1,ag,bg));

lamm=0; mn=1; TotExpProfitRate=zeros(length(p1),length(p2),length(cc));


for c_incr=1:length(cc) %
    lamm=lamm+1;
    
    for jj=1:length(p2)
        for ii=1:length(p1) %
            
            if lambda1(p1(ii),p2(jj),lambda) <mu1(p1(ii),p2(jj)) % check transp. queue stability condition
                ER = cc(c_incr)*(1-2*gammainc(lambda2(p1(ii),p2(jj),lambda)/mu2,theta,'upper'));
%                 ER = cc*(1-2*(igamma(theta(theta_incr),lambda*(p2(jj)/(g1*p1(ii)+p2(jj)))))/gamma(theta(theta_incr)));
                TotExpProfitRate(ii,jj,lamm)=((1-g1)*lambda1(p1(ii),p2(jj),lambda)*p1(ii)+ ER -lambda2(p1(ii),p2(jj),lambda)*p2(jj))/mn;
                
            else % if transp. queue stability condition is not met, announce nothing
                TotExpProfitRate(ii,jj,lamm)=NaN;
                
            end
        end
    end
end

figure
for LL=1:size(TotExpProfitRate,3);
    surf(p1,p2,TotExpProfitRate(:,:,10)');
    hold on
end
xlabel('p_1')
ylabel('p_2')
zlabel('revenue rate')
title('Total expected profit rate vs c = a representative c')


% code to plot the optimal prices versus theta
LLL=[]; LINEARINDEX=[]; P1=[]; P2=[];
% here it extracts the optimal prices from TotExpProfitRate
for LL=1:size(TotExpProfitRate,3);
    [maxValue, linearIndexesOfMaxes] = max(max(TotExpProfitRate(:,:,LL)));
    if min(min(isnan((TotExpProfitRate(:,:,LL)))))~=1 % only the non nan exp. rates (they have at least one max) are kept
        [rowsOfMaxes, colsOfMaxes] = find(TotExpProfitRate(:,:,LL) == maxValue);
        LINEARINDEX=[LINEARINDEX; sub2ind(size(TotExpProfitRate), rowsOfMaxes, colsOfMaxes, LL)];
        P1=[P1; p1(rowsOfMaxes)];
        P2=[P2; p2(colsOfMaxes)];
        LLL=[LLL; LL];
    end
end

%%% plot - maximizing prices vs theta
figure
scatter(cc(LLL)/mn,P1,[],TotExpProfitRate(LINEARINDEX),'filled','LineWidth',1.8)
hold on
scatter(cc(LLL)/mn,P2,[],TotExpProfitRate(LINEARINDEX),'o','LineWidth',1.8,'MarkerEdgeColor',[0 .5 .5]) % ,
xlabel('c')
ylabel('prices')
legend('p1*','p2*')
title('Prices for max profit vs c')
grid on

%%% plot - Total profit rate at max prices vs lambda
figure
plot(cc(LLL),TotExpProfitRate(LINEARINDEX),'-o','LineWidth',2)
xlabel('c')
ylabel('revenue rate')
title('Total expected profit rate @ max prices vs c')