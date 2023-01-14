close all; clear all; clc

set(0,'DefaultFigureUnits','Inches');
set(0,'DefaultAxesUnits','Inches');

%% Import Data

rootDir = dir('G:\Pressure\rampJuly24\averagedData');
set = {'(1)(2)(5)(10)'; '(3)(6)(9)(11)(14)'; '(4)(7)(8)(12)(13)'};

A20S1 = cell(1,3);
A20S2 = cell(1,3);
A30S1 = cell(1,3);
A30S2 = cell(1,3);

for file = 1:3

      A20S1{file}  =  xlsread(['averagedPressure_', ... 
          sprintf('%s',set{file}),'_20deg_sStar1'], 'averagedValues');
      A20S2{file}  =  xlsread(['averagedPressure_', ... 
          sprintf('%s',set{file}),'_20deg_sStar2'], 'averagedValues');
      A30S1{file}  =  xlsread(['averagedPressure_', ... 
          sprintf('%s',set{file}),'_30deg_sStar1'], 'averagedValues'); 
      A30S2{file}  =  xlsread(['averagedPressure_', ... 
          sprintf('%s',set{file}),'_30deg_sStar2'], 'averagedValues');
      
      lgth1(file) = length(A20S1{file});
      lgth2(file) = length(A20S2{file});
      lgth3(file) = length(A30S1{file});
      lgth4(file) = length(A30S2{file});
      
end

for set = 1:4
D(set) = min(sprintf('lgth%d',set));
end

L = min(D);

pSet1 = [A20S1{1}(1:D) A20S1{2}(1:D) A20S1{2}(1:D)];

%% Experimental Parameters

sampleRate = 1000;
NumP  = 5;
rho = 998;
c=0.3;
U=1.0; %velocity
T_s=1/sampleRate; %sampling period (s)
vis=1.004e-6; %kinematic viscosity of water at 20deg C (m^2/s)
Re=U*c/vis; %Reynolds number

%% Non-dimensional time and distance
t =    {linspace(0,T_s*length(deg10(:,1)),length(deg10(:,1))), ...
       linspace(0,T_s*length(deg20(:,1)),length(deg20(:,1))), ...
       linspace(0,T_s*length(deg30(:,1)),length(deg30(:,1)))};
 
t_star = {t{1}*U/c t{2}*U/c t{3}*U/c};
s_star = {t_star{1}/2 t_star{2}/2 t_star{3}/2};

%% Smooth Values and Calulate Cp

P_Range = 17000; %pressure in Pa
V_Range = 5;


deg10_P = zeros(length(deg10(:,1)),NumP);
deg20_P = zeros(length(deg20(:,1)),NumP);
deg30_P = zeros(length(deg30(:,1)),NumP);

Cp10 = zeros(length(deg10(:,1)),NumP);
Cp20 = zeros(length(deg20(:,1)),NumP);
Cp30 = zeros(length(deg30(:,1)),NumP);

for i = 1:NumP
    deg10(:,i) = smoother200(deg10(:,i));
    deg10_P(:,i) = deg10(:,i).*P_Range./V_Range;
    deg20(:,i) = smoother200(deg20(:,i));
    deg20_P(:,i) = deg20(:,i).*P_Range./V_Range;
    deg30(:,i) = smoother200(deg30(:,i));
    deg30_P(:,i) = deg30(:,i).*P_Range./V_Range;

    dH10 = mean(deg10_P(1:500,i));
    dH20 = mean(deg20_P(1:500,i));
    dH30 = mean(deg30_P(1:500,i));

    if i == 1
        
        Cp10(:,i) = -1*(deg10_P(:,i) - dH10)/(0.5*rho*U^2);
        Cp20(:,i) = -1*(deg20_P(:,i) - dH20)/(0.5*rho*U^2);
        Cp30(:,i) = -1*(deg30_P(:,i) - dH30)/(0.5*rho*U^2);
        
    elseif i == 2 || i == 3
        
        Cp10(:,i) = Cp10(:,i-1) + (deg10_P(:,i) - dH10)/(0.5*rho*U^2);
        Cp20(:,i) = Cp20(:,i-1) + (deg20_P(:,i) - dH20)/(0.5*rho*U^2);
        Cp30(:,i) = Cp30(:,i-1) + (deg30_P(:,i) - dH30)/(0.5*rho*U^2);
        
    else
        
        Cp10(:,i) = Cp10(:,i-1) - (deg10_P(:,i) - dH10)/(0.5*rho*U^2);
        Cp20(:,i) = Cp20(:,i-1) - (deg20_P(:,i) - dH20)/(0.5*rho*U^2);
        Cp30(:,i) = Cp30(:,i-1) - (deg30_P(:,i) - dH30)/(0.5*rho*U^2);
        
    end

end

%% Calculate mean pressure over the range sLow < s < sHigh

sLow = 8;
sHigh = 18;

lb =  [find(s_star{1} > sLow, 1) ...
       find(s_star{2} > sLow, 1) ...
       find(s_star{3} > sLow, 1) ...
       ];
   

ub =  [find(s_star{1} > sHigh, 1) ...
       find(s_star{2} > sHigh, 1) ...
       find(s_star{3} > sHigh, 1) ...
       ];

Cp10_mean = zeros(1,NumP);
Cp20_mean = zeros(1,NumP);
Cp30_mean = zeros(1,NumP);

for j = 1:NumP
    Cp10_mean(j) = mean(Cp10(lb(1):ub(1),j));
    Cp20_mean(j) = mean(Cp20(lb(1):ub(1),j));
    Cp30_mean(j) = mean(Cp30(lb(1):ub(1),j));
end

%% Plot Steady Pressure Distribution

fig=figure(2);
set(fig,'position',[1.5 1.5 10 6],'color',[1 1 1])
set(gca,'fontsize',13)

hold on
% plot(x_c, -1*Cp_mean10, '-.oc', x_c, -1*Cp_mean20, '-.xm', x_c, ...
%     -1*Cp_mean30, '-.^g',xFoil(3:end,1),xFoil(3:end,2),'-k', ... 
%     'markersize',10,'linewidth',2)
plot(x_c, -1*Cp10_mean, '-.oc', x_c, -1*Cp20_mean, '-.xm', x_c, ...
    -1*Cp30_mean, '-.^g','markersize',10,'linewidth',2)

xlabel('$x/c$','Interpreter','Latex','FontSize', 17);
ylabel('$-\langle C_p \rangle$','Interpreter','Latex','FontSize',17);
% h = legend('$10^{\circ}$','$20^{\circ}$','$30^{\circ}$','NACA 0012 2D Airfoil');
h = legend('$10^{\circ}$','$20^{\circ}$','$30^{\circ}$');
set(h,'Interpreter','Latex','FontSize',13)

% axis([0 1 -0.3 2])

hold on
grid on
ax = gca;
ax.GridLineStyle = '--';
