close all
clear all
clc

set(0,'DefaultFigureUnits','Inches');
set(0,'DefaultAxesUnits','Inches');

%%%%% Disclaimer - this code is pure unadulterated monkey shit, please take
%%%%% inspiration from my code as how not too do things. The VBA code is
%%%%% fast and works really well; if your useing this, just delete
%%%%% everything after the data import.


%% Import Data
filename10 = strcat('C:\Users\Matt Marzanek\Desktop\pressureMeasurements\', ...
    '\pressureMeasurements\steadyJuly11\leadingEdge\10deg\', ...
    'averagedPressure_(1)(4)(9)(16)_10deg_steady.xlsx');

filename20 = strcat('C:\Users\Matt Marzanek\Desktop\pressureMeasurements\', ...
    '\pressureMeasurements\steadyJuly11\leadingEdge\20deg\', ...
    'averagedPressure_(1)(4)(9)(16)_20deg_steady.xlsx');
    
filename30 = strcat('C:\Users\Matt Marzanek\Desktop\pressureMeasurements\', ...
    '\pressureMeasurements\steadyJuly11\leadingEdge\30deg\', ...
    'averagedPressure_(1)(4)(9)(16)_30deg_steady.xlsx');


deg10 = xlsread(filename10,'averagedValues');
deg20 = xlsread(filename20,'averagedValues');
deg30 = xlsread(filename30,'averagedValues');


%% Experimental Parameters

sampleRate = 1000;
refHead = 4000;
rho = 998;
c=0.3;
area=0.09;
U=1.0; %velocity
D=11; %distance travelled by carriage (m)
T_s=1/sampleRate; %sampling period (s)
vis=1.004e-6; %kinematic viscosity of water at 20deg C (m^2/s)

t_accel = 6; %t_star
s_accel=t_accel./2; %s_star
Re=U*c/vis; %corresponding Reynolds number
a=(1./(t_accel*c)).*(U^2); %vector of accelerations

%% Calculate dimentionless distances

accel_d = U^2/(2*a); %distance for accel
accel_t = U/a;
steady_d = D-(2*accel_d); %remaining distance
steady_t = steady_d/U;

t10 = linspace(0,T_s*length(deg10(:,1)),length(deg10(:,1)));
t20 = linspace(0,T_s*length(deg20(:,1)),length(deg20(:,1)));
t30 = linspace(0,T_s*length(deg30(:,1)),length(deg30(:,1)));

t_star10 = t10*U/c;
t_star20 = t20*U/c;
t_star30 = t30*U/c;

s_star10 = t_star10/2;
s_star20 = t_star20/2;
s_star30 = t_star30/2;

%% Calculate Pressure and Smooth Values

P_Range = 17000; %pressure in Pa
V_Range = 5;

for i = 1:length(deg10(1,:))
deg10(:,i) = smoother200(deg10(:,i));
deg10_P(:,i) = deg10(:,i).*P_Range./V_Range;
deg20(:,i) = smoother200(deg20(:,i));
deg20_P(:,i) = deg20(:,i).*P_Range./V_Range;
deg30(:,i) = smoother200(deg30(:,i));
deg30_P(:,i) = deg30(:,i).*P_Range./V_Range;
end

% should have used a loop for all of this :( 

%dH = [0.2*c*sind(10), 0.2*c*sind(20), 0.2*c*sind(30)].*(rho*9.81);

Cp110 = -1*(deg10_P(:,1) - mean(deg10_P(1:10,1)))/(0.5*rho*U^2);
Cp120 = -1*(deg20_P(:,1) - mean(deg20_P(1:10,1)))/(0.5*rho*U^2);
Cp130 = -1*(deg30_P(:,1) - mean(deg30_P(1:10,1)))/(0.5*rho*U^2);

Cp210 = Cp110 + (deg10_P(:,2) - mean(deg10_P(1:10,2)))/(0.5*rho*U^2);
Cp220 = Cp120 + (deg20_P(:,2) - mean(deg20_P(1:10,2)))/(0.5*rho*U^2);
Cp230 = Cp130 + (deg30_P(:,2) - mean(deg30_P(1:10,2)))/(0.5*rho*U^2);

Cp310 = Cp210 + (deg10_P(:,3) - mean(deg10_P(1:10,3)))/(0.5*rho*U^2);
Cp320 = Cp220 + (deg20_P(:,3) - mean(deg20_P(1:10,3)))/(0.5*rho*U^2);
Cp330 = Cp230 + (deg30_P(:,3) - mean(deg30_P(1:10,3)))/(0.5*rho*U^2);

Cp410 = Cp310 - (deg10_P(:,4) - mean(deg10_P(1:10,4)))/(0.5*rho*U^2);
Cp420 = Cp320 - (deg20_P(:,4) - mean(deg20_P(1:10,4)))/(0.5*rho*U^2);
Cp430 = Cp330 - (deg30_P(:,4) - mean(deg30_P(1:10,4)))/(0.5*rho*U^2);

Cp510 = Cp410 - (deg10_P(:,5) - mean(deg10_P(1:10,5)))/(0.5*rho*U^2);
Cp520 = Cp420 - (deg20_P(:,5) - mean(deg20_P(1:10,5)))/(0.5*rho*U^2);
Cp530 = Cp430 - (deg30_P(:,5) - mean(deg30_P(1:10,5)))/(0.5*rho*U^2);

%% Steady State

lb10 = find(s_star10>8,1);
ub10 = find(s_star10>18,1);

lb20 = find(s_star20>8,1);
ub20 = find(s_star20>18,1);

lb30 = find(s_star30>8,1);
ub30 = find(s_star30>18,1);

Cp110_mean = mean(Cp110(lb10:ub10));
Cp210_mean = mean(Cp210(lb10:ub10));
Cp310_mean = mean(Cp310(lb10:ub10));
Cp410_mean = mean(Cp410(lb10:ub10));
Cp510_mean = mean(Cp510(lb10:ub10));
Cp_mean10 = [Cp110_mean Cp210_mean Cp310_mean Cp410_mean Cp510_mean];

Cp120_mean = mean(Cp120(lb20:ub20));
Cp220_mean = mean(Cp220(lb20:ub20));
Cp320_mean = mean(Cp320(lb20:ub20));
Cp420_mean = mean(Cp420(lb20:ub20));
Cp520_mean = mean(Cp520(lb10:ub10));
Cp_mean20 = [Cp120_mean Cp220_mean Cp320_mean Cp420_mean Cp520_mean];

Cp130_mean = mean(Cp130(lb30:ub30));
Cp230_mean = mean(Cp230(lb30:ub30));
Cp330_mean = mean(Cp330(lb30:ub30));
Cp430_mean = mean(Cp430(lb30:ub30));
Cp530_mean = mean(Cp530(lb10:ub10));
Cp_mean30 = [Cp130_mean Cp230_mean Cp330_mean Cp430_mean Cp530_mean];

x_c = [0.1 0.2 0.3 0.4 0.5];

% RMS
% Cp10sqr = [(Cp110_mean-Cp110(lb10:ub10)).^2, ...
%            (Cp210_mean-Cp210(lb10:ub10)).^2, ...
%            (Cp310_mean-Cp310(lb10:ub10)).^2,  ...
%            (Cp410_mean-Cp410(lb10:ub10)).^2];
%        
% Cp20sqr = [(Cp120_mean-Cp120(lb20:ub20)).^2, ...
%            (Cp220_mean-Cp220(lb20:ub20)).^2, ...
%            (Cp320_mean-Cp320(lb20:ub20)).^2,  ...
%            (Cp420_mean-Cp420(lb20:ub20)).^2];
%        
% Cp30sqr = [(Cp130_mean-Cp130(lb30:ub30)).^2, ...
%            (Cp230_mean-Cp230(lb30:ub30)).^2, ...
%            (Cp330_mean-Cp330(lb30:ub30)).^2,  ...
%            (Cp430_mean-Cp430(lb30:ub30)).^2];

%peak to peak
Cp10sqr = [range(Cp110(lb10:ub10)), ...
           range(Cp210(lb10:ub10)), ...
           range(Cp310(lb10:ub10)),  ...
           range(Cp410(lb10:ub10))];
       
Cp20sqr = [range(Cp120(lb20:ub20)), ...
           range(Cp220(lb20:ub20)), ...
           range(Cp320(lb20:ub20)),  ...
           range(Cp420(lb20:ub20))];
       
Cp30sqr = [range(Cp130(lb30:ub30)), ...
           range(Cp230(lb30:ub30)), ...
           range(Cp330(lb30:ub30)),  ...
           range(Cp430(lb30:ub30))];

for i = 1:4
%     Cp10RMS(i) = sqrt(sum(Cp10sqr(i,:))/length(lb10:ub10));
%     Cp20RMS(i) = sqrt(sum(Cp20sqr(i,:))/length(lb20:ub20));
%     Cp30RMS(i) = sqrt(sum(Cp30sqr(i,:))/length(lb30:ub30));
    Cp10RMS(i) = Cp10sqr(i)/2;
    Cp20RMS(i) = Cp20sqr(i)/2;
    Cp30RMS(i) = Cp30sqr(i)/2;
end

%% Calculate FFT

%[FFT] = calculateFFT(sampleRate, lb10, ub10, lb20, ub20, lb30, ub30);
R10 = length(lb10:ub10);
R20 = length(lb20:ub20);
R30 = length(lb30:ub30);

f10 = sampleRate*(0:(R10/2))./R10;
f20 = sampleRate*(0:(R20/2))./R20;
f30 = sampleRate*(0:(R30/2))./R30;

dtnd = detrend(Cp110(lb10:ub10));
FFT110 = fft(dtnd,R10);
temp = FFT110.*conj(FFT110)/R10;
FFT110 = temp(1:R10/2+1);
FFT110(2:end-1) = 2*FFT110(2:end-1);

dtnd = detrend(Cp210(lb10:ub10));
FFT210 = fft(dtnd,R10);
temp = FFT210.*conj(FFT210)/R10;
FFT210 = temp(1:R10/2+1);
FFT210(2:end-1) = 2*FFT210(2:end-1);

dtnd = detrend(Cp310(lb10:ub10));
FFT310 = fft(dtnd,R10);
temp = FFT310.*conj(FFT310)/R10;
FFT310 = temp(1:R10/2+1);
FFT310(2:end-1) = 2*FFT310(2:end-1);

dtnd = detrend(Cp410(lb10:ub10));
FFT410 = fft(dtnd,R10);
temp = FFT410.*conj(FFT410)/R10;
FFT410 = temp(1:R10/2+1);
FFT410(2:end-1) = 2*FFT410(2:end-1);

dtnd = detrend(Cp120(lb20:ub20));
FFT120 = fft(dtnd,R20);
temp = FFT120.*conj(FFT120)/R20;
FFT120 = temp(1:R20/2+1);
FFT120(2:end-1) = 2*FFT120(2:end-1);

dtnd = detrend(Cp220(lb20:ub20));
FFT220 = fft(dtnd,R20);
temp = FFT220.*conj(FFT220)/R20;
FFT220 = temp(1:R20/2+1);
FFT220(2:end-1) = 2*FFT220(2:end-1);

dtnd = detrend(Cp320(lb20:ub20));
FFT320 = fft(dtnd,R20);
temp = FFT320.*conj(FFT320)/R20;
FFT320 = temp(1:R20/2+1);
FFT320(2:end-1) = 2*FFT320(2:end-1);

dtnd = detrend(Cp420(lb20:ub20));
FFT420 = fft(dtnd,R20);
temp = FFT420.*conj(FFT420)/R20;
FFT420 = temp(1:R20/2+1);
FFT420(2:end-1) = 2*FFT420(2:end-1);

dtnd = detrend(Cp130(lb30:ub30));
FFT130 = fft(dtnd,R20);
temp = FFT130.*conj(FFT130)/R30;
FFT130 = temp(1:R30/2+1);
FFT130(2:end-1) = 2*FFT130(2:end-1);

dtnd = detrend(Cp230(lb30:ub30));
FFT230 = fft(dtnd,R20);
temp = FFT230.*conj(FFT230)/R30;
FFT230 = temp(1:R30/2+1);
FFT230(2:end-1) = 2*FFT230(2:end-1);

dtnd = detrend(Cp330(lb30:ub30));
FFT330 = fft(dtnd,R20);
temp = FFT330.*conj(FFT330)/R30;
FFT330 = temp(1:R30/2+1);
FFT330(2:end-1) = 2*FFT330(2:end-1);

dtnd = detrend(Cp430(lb30:ub30));
FFT430 = fft(dtnd,R20);
temp = FFT430.*conj(FFT430)/R30;
FFT430 = temp(1:R30/2+1);
FFT430(2:end-1) = 2*FFT430(2:end-1);

FFT = [FFT110 FFT210 FFT310 FFT410 FFT120 FFT220 FFT320 FFT420 ...
        FFT130 FFT230 FFT330 FFT430];

%% Plot

fig=figure(1);
set(fig,'position',[1 1 10 6],'color',[1 1 1])
set(gca,'fontsize',13)
% ax = gca;
% ax.GridLineStyle = '--';
% ax.LineWidth = 2;
%set(gca,'position',[0.65 0.6 6.4 4.6])

hold on
grid on

subplot(2,2,1)
plot(s_star10, Cp110, 'c', ...
     s_star20, Cp120, 'm', ...
     s_star30, Cp130, 'g','linewidth',1.35);
xlabel('$s^*$','Interpreter','Latex','FontSize', 17);
ylabel('$C_p$','Interpreter','Latex','FontSize',17);
axis([0 20 -4 1])
grid on
ax = gca;
ax.GridLineStyle = '--';

subplot(2,2,2)
plot(s_star10, Cp210, 'c', ...
     s_star20, Cp220, 'm', ...
     s_star30, Cp330, 'g','linewidth',1.35);
xlabel('$s^*$','Interpreter','Latex','FontSize', 17);
ylabel('$C_p$','Interpreter','Latex','FontSize',17);
h = legend('$10^{\circ}$','$20^{\circ}$','$30^{\circ}$');
set(h,'Interpreter','Latex','FontSize',13)
axis([0 20 -1.5 0.5])
grid on
ax = gca;
ax.GridLineStyle = '--';

subplot(2,2,3)
plot(s_star10, Cp310, 'c', ...
     s_star20, Cp320, 'm', ...
     s_star30, Cp330, 'g','linewidth',1.35);
xlabel('$s^*$','Interpreter','Latex','FontSize', 17);
ylabel('$C_p$','Interpreter','Latex','FontSize',17);
axis([0 20 -1.5 0.5])
grid on
ax = gca;
ax.GridLineStyle = '--';

subplot(2,2,4)
plot(s_star10, Cp410, 'c', ...
     s_star20, Cp420, 'm', ...
     s_star30, Cp430, 'g','linewidth',1.35);
xlabel('$s^*$','Interpreter','Latex','FontSize', 17);
ylabel('$C_p$','Interpreter','Latex','FontSize',17);
axis([0 20 -1.5 0.5])
grid on
ax = gca;
ax.GridLineStyle = '--';

%**********************Steady************************%

xFoil = xlsread('nacaXfoil.csv');
fig=figure(2);
set(fig,'position',[1.5 1.5 10 6],'color',[1 1 1])
set(gca,'fontsize',13)
hold on
plot(x_c, -1*Cp_mean10, '-.oc', x_c, -1*Cp_mean20, '-.xm', x_c, ...
    -1*Cp_mean30, '-.^g',xFoil(3:end,1),xFoil(3:end,2),'-k', ... 
    'markersize',10,'linewidth',2)
plot(x_c, -1*Cp_mean10, '-.oc', x_c, -1*Cp_mean20, '-.xm', x_c, ...
    -1*Cp_mean30, '-.^g','markersize',10,'linewidth',2)
xlabel('$x/c$','Interpreter','Latex','FontSize', 17);
ylabel('$-\langle C_p \rangle$','Interpreter','Latex','FontSize',17);
h = legend('$10^{\circ}$','$20^{\circ}$','$30^{\circ}$','NACA 0012 2D Airfoil');
set(h,'Interpreter','Latex','FontSize',13)

axis([0 1 -0.3 2])

hold on
grid on
ax = gca;
ax.GridLineStyle = '--';

%%  FFT Plots

fig = figure(3);
set(fig,'position',[1.2 1.2 10 6],'color',[1 1 1])
set(gca,'fontsize',13)


subplot(2,2,1)
semilogx(f30, FFT130, 'k')
xlabel('$f$ (Hz)','Interpreter','Latex','FontSize', 17);
ylabel('$P.S.D$','Interpreter','Latex','FontSize',17);

%axis([0 100 0 0.7])

subplot(2,2,2)
semilogx(f30, FFT230, 'k')
xlabel('$f$ (Hz)','Interpreter','Latex','FontSize', 17);
ylabel('$P.S.D$','Interpreter','Latex','FontSize',17);


%axis([0 100 0 0.7])

subplot(2,2,3)
semilogx(f30, FFT330, 'k')
xlabel('$f$ (Hz)','Interpreter','Latex','FontSize', 17);
ylabel('$P.S.D$','Interpreter','Latex','FontSize',17);


%axis([0 100 0 0.5])

subplot(2,2,4)
semilogx(f30, FFT430, 'k')
xlabel('$f$ (Hz)','Interpreter','Latex','FontSize', 17);
ylabel('$P.S.D$','Interpreter','Latex','FontSize',17);


%axis([0 100 0 0.5])

%% RAW DATA

fig=figure(4);
set(fig,'position',[1 1 10 6],'color',[1 1 1])
set(gca,'fontsize',13)
% ax = gca;
% ax.GridLineStyle = '--';
% ax.LineWidth = 2;
%set(gca,'position',[0.65 0.6 6.4 4.6])

hold on
grid on

subplot(2,2,1)
plot(s_star10, deg10_P(:,1), 'c', ...
     s_star20, deg20_P(:,1), 'm', ...
     s_star30, deg30_P(:,1), 'g','linewidth',1.35);
xlabel('$s^*$','Interpreter','Latex','FontSize', 17);
ylabel('$C_p$','Interpreter','Latex','FontSize',17);
%axis([0 20 -3 0])
grid on
ax = gca;
ax.GridLineStyle = '--';

subplot(2,2,2)
plot(s_star10, deg10_P(:,2), 'c', ...
     s_star20, deg20_P(:,2), 'm', ...
     s_star30, deg30_P(:,2), 'g','linewidth',1.35);
xlabel('$s^*$','Interpreter','Latex','FontSize', 17);
ylabel('$C_p$','Interpreter','Latex','FontSize',17);
h = legend('$10^{\circ}$','$20^{\circ}$','$30^{\circ}$');
set(h,'Interpreter','Latex','FontSize',13)
%axis([0 20 -0.5 0.5])
grid on
ax = gca;
ax.GridLineStyle = '--';

subplot(2,2,3)
plot(s_star10, deg10_P(:,3), 'c', ...
     s_star20, deg20_P(:,3), 'm', ...
     s_star30, deg30_P(:,3), 'g','linewidth',1.35);
xlabel('$s^*$','Interpreter','Latex','FontSize', 17);
ylabel('$C_p$','Interpreter','Latex','FontSize',17);
%axis([0 20 -0.5 0.5])
grid on
ax = gca;
ax.GridLineStyle = '--';

subplot(2,2,4)
plot(s_star10, deg10_P(:,4), 'c', ...
     s_star20, deg20_P(:,4), 'm', ...
     s_star30, deg30_P(:,4), 'g','linewidth',1.35);
xlabel('$s^*$','Interpreter','Latex','FontSize', 17);
ylabel('$C_p$','Interpreter','Latex','FontSize',17);
%axis([0 20 -0.5 0.5])
grid on
ax = gca;
ax.GridLineStyle = '--';
