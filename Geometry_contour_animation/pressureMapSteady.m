close all; clear all; clc;
set(0,'DefaultFigureUnits','Inches');
set(0,'DefaultAxesUnits','Inches');
%set(0, 'DefaultFigurePosition', [14.4375   -2.3958  10.0000    6.0000])

%% Load Data

load('alpha20_sStar1.mat')
load('alpha20_sStar2.mat')
load('alpha30_sStar1.mat')
load('alpha30_sStar2.mat')
load('Length.mat')

%% Input Variables
U = 1;
Uf = 1.5;
c = 0.3;
sampleRate = 1000;
T_s = 1/sampleRate;
t = linspace(0,T_s*L,L);
t_star = t*U/c;
%s_star = t_star/2;

t_star = t_star - t_star(7600);
t_star(1:7600) = 0;

%% Trill shit dawg.. seriously don't bother with this crap here
sRamp = [1,2];
a = U*(Uf-U)./(sRamp.*c) + (Uf-U)^2./(2*sRamp*c);
t3 = (Uf - U)./a;
t4 = [4 4]/Uf;%dSSR/Uf;
for i = 1:2
    SR{i,:} = linearKinematics(U,T_s,a(i),sRamp(i),c,t3(i));
    SRS{i,:} = linspace(max(SR{i,:}),max(SR{i,:})+4/c,t4(i)/T_s);
end
    SStr = {[SR{1,:}' SRS{1,:}] [SR{2,:}' SRS{2,:}]};

%% Movie

%  writerObj = VideoWriter('pressureAndForces30deg_sStar2.avi');
%  open(writerObj);

%% Generate Pressure Map

x = [0.1 0.1 0.1 0.1 0.2 0.2 0.5 0.2 0.5 0.3 0.3 0.4 0.3 0.4];  
y = [0.2 0.4 0.6 0.8 0.4 0.6 0.6 0.8 0.8 0.4 0.6 0.6 0.8 0.8];

for i = 8200%5501:6:10080 %1:12:10080 %5500:10:10300 %7840 8080
    
tStep1 = 5500; %find(s_star > 10,1);
tStep2 = 7600;%find(s_star > 18,1);

p1 = mean(Cp1(tStep1:tStep2,:));
p2 = mean(Cp2(tStep1:tStep2,:));
p3 = mean(Cp3(tStep1:tStep2,:));
p4 = mean(Cp4(tStep1:tStep2,:));

% if i > 7840
% d = SStr{1,1}(1,i-7840) - SStr{1,1}(1,i-7839);
% p1 = Cp1(i,:)./((d*c/T_s)^2);
% p2 = Cp2(i,:)./((d*c/T_s)^2);
% p3 = Cp3(i,:)./((d*c/T_s)^2);
% p4 = Cp4(i,:)./((d*c/T_s)^2);
% else
% p1 = Cp1(i,:);
% p2 = Cp2(i,:);
% p3 = Cp3(i,:);
% p4 = Cp4(i,:);
% end

p1 = [p1(:,1:4),p1(:,6:end)];
p2 = [p2(:,1:4),p2(:,6:end)];
p3 = [p3(:,1:4),p3(:,6:end)];
p4 = [p4(:,1:4),p4(:,6:end)];

% for i = 1:length(Cp1(1,:))
% p1(i) = pSample(i);
% end
             
%[xG, yG] = meshgrid(-0.01:0.001:1, -0.01:0.001:1);
chord = 1;
grSize = 1200;
[xx,yy]=meshgrid(linspace(-1,1,grSize),linspace(0,1,grSize));

vG(:,:,1) = griddata(x,y,p1,xx,yy,'cubic');
vG(:,:,2) = griddata(x,y,p2,xx,yy,'cubic');
vG(:,:,3) = griddata(x,y,p3,xx,yy,'cubic');
vG(:,:,4) = griddata(x,y,p4,xx,yy,'cubic');

%% Colormap
l=50;
n=50;
m=75;

segment1=...
    [   linspace(0/255,1,l)
        linspace(0/255,1,l)
        linspace(255/255,1,l)];
    
segment2=...
    [   linspace(1,1,m)
        linspace(1,1,m)
        linspace(1,1,m)];

segment3=...
    [   linspace(1,204/255,l)
        linspace(1,0,l)
        linspace(1,0,l)];
        
shittyBlue = [segment1,segment2]';
shittyBlueTwo = [segment1,segment2,segment3]'; %electric bugaloo

load('colorbar.mat')

%% Plot pressure map

P = [p1; p2; p3; p4];

plotWingSteady(xx,yy,vG,grSize,3,[-1.5 0],i,SStr)
%plotPressureMaps(xG, yG, vG, x, y, P, [-1.5 0],i,t_star,SStr)
% magma=magma();
% cmap = magma;%four_color_vort;
colormap(shittyBlue)
%colormap(cmap(1:112,:)); %only uses the cold portion of the color map
% colormap(cmap(112:225,:)); %only uses the warmportion of the color map
colormap jet
colorbar off

% %cd('D:\Critical_Data_Backup\1_Forces\rampForces\Plotting')
% 
% axes('Position',[5.65 3 4 2.5])
% box on
% shiftRamp
% cd('D:\Critical_Data_Backup\2_Pressure\rampJuly24\averagedData')

% if i < 7840
%     t = annotation('textbox');
%     t.String = 'Steady State';
%     t.Interpreter = 'latex';
%     t.FontSize = 32;
%     t.LineStyle = 'none';
%     t.Position = [0.2 0.75 0.1 0.1];
%     else
%     str = num2str(SStr{1,1}(1,i-7840),'%0.2f');
%     t = annotation('textbox');
%     t.String = ['$s^*$ =', str];
%     t.Interpreter = 'latex';
%     t.FontSize = 32;
%     t.LineStyle = 'none';
%     t.Position = [0.2 0.75 0.1 0.1];
%     yS = linspace(-0.8,3,100);
%     s1 = ones(length(yS),1)*SStr{1,1}(1,i-7840);
%     plot(s1,yS,':k','linewidth',3)
% 
%     end
    
    %% Make Movie
     %N = myaa;
%      temp =getframe(gcf);
%      writeVideo(writerObj,temp);
%     disp(i)
%     close all
%     
end

   % close(writerObj);

%% Plot pressure traces
% t_star = t*U/c;
% 
% fig=figure(2);
% set(fig,'position',[1 1 8 8],'color',[1 1 1])
% set(gca,'fontsize',13)
% 
% ax = gca;
% ax.GridLineStyle = '--';
% ax.LineWidth = 2;
% set(gca,'position',[0.65 0.6 6.4 4.6])
% 
% hold on
% grid on
% s_star = s_star';
% 
% subplot(2,2,1)
% plot(t_star, Cp1(:,5), 'c', ...
%      t_star, Cp2(:,5), 'm', ...
%      t_star, Cp3(:,5), 'g', ...
%      t_star, Cp4(:,5), 'r',  'linewidth',1.35);
%  
% xlabel('$s^*$','Interpreter','Latex','FontSize', 17);
% ylabel('$C_p$','Interpreter','Latex','FontSize',17);
% axis([0 40 -4 1])
% grid on
% ax = gca;
% ax.GridLineStyle = '--';
% 
% subplot(2,2,2)
% plot(t_star, Cp1(:,10), 'c', ...
%      t_star, Cp2(:,10), 'm', ...
%      t_star, Cp3(:,10), 'g', ...
%      t_star, Cp4(:,10), 'r',  'linewidth',1.35);
%  
% xlabel('$s^*$','Interpreter','Latex','FontSize', 17);
% ylabel('$C_p$','Interpreter','Latex','FontSize',17);
% h = legend('$10^{\circ}$','$20^{\circ}$','$30^{\circ}$');
% set(h,'Interpreter','Latex','FontSize',13)
% axis([0 40 -4 1])
% grid on
% ax = gca;
% ax.GridLineStyle = '--';
% 
% subplot(2,2,3)
% plot(t_star, Cp1(:,12), 'c', ...
%      t_star, Cp2(:,12), 'm', ...
%      t_star, Cp3(:,12), 'g', ...
%      t_star, Cp4(:,12), 'r',  'linewidth',1.35);
%  
% xlabel('$s^*$','Interpreter','Latex','FontSize', 17);
% ylabel('$C_p$','Interpreter','Latex','FontSize',17);
% axis([0 40 -4 1])
% grid on
% ax = gca;
% ax.GridLineStyle = '--';
% 
% subplot(2,2,4)
% plot(t_star, Cp1(:,13), 'c', ...
%      t_star, Cp2(:,13), 'm', ...
%      t_star, Cp3(:,13), 'g', ...
%      t_star, Cp4(:,13), 'r',  'linewidth',1.35);
% xlabel('$s^*$','Interpreter','Latex','FontSize', 17);
% ylabel('$C_p$','Interpreter','Latex','FontSize',17);
% axis([0 40 -4 1])
% grid on
% ax = gca;
% ax.GridLineStyle = '--';

% D = 11.5;
% si = 3;
% ai = 5/9;
% sf = 6.75;
% af = -5/9;
% sRamp = [1,2];
% a = U*(Uf-U)./(sRamp.*c) + (Uf-U)^2./(2*sRamp*c);
% 
%     dSS = D - (si + sf + sRamp)*c;
%     dSSI = (6-si*c);
%     dSSR = dSS - dSSI;
%     t1 =  U/ai;
%     t2 = dSSI/U;
%     t3 = (Uf - U)./a;
%     t4 = [4 2]*c/Uf;%dSSR/Uf;
%     t5 = abs(Uf/af);
%     tEst = t1 + t2 + t3 +  t4 + t5;
% 
% s1 = linearKinematics(0,T_s,ai,si,c,t1);
% s2 = linspace(max(s1),dSSI/c+max(s1),round(t2/T_s));
% s5 = linearKinematics(Uf,T_s,af,sf,c,t5);
% 
% for i = 1:2
%     s3 = linearKinematics(U,T_s,a(i),sRamp(i),c,t3(i));
%     sabs3 = s3 + max(s2);
%     s4 = linspace(max(sabs3), dSSR(i)/c+max(sabs3),t4(i)/T_s);
%     sabs5 = s5 + max(s4);
%     s{i,:} = [s1;s2';sabs3;s4';sabs5]';
%     Lth(i) = length(s{i,:});
% end
