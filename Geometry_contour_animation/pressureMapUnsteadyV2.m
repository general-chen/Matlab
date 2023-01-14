close all; clear all; clc;
set(0,'DefaultFigureUnits','Inches');
set(0,'DefaultAxesUnits','Inches');

%%%************************************************************************
%**************************************************************************

%V2: moved plotting to separate function to reduce clutter, now plots as 4
%separate sub-plots with annotation and titles
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
s_star = t_star/2;

t_star = t_star - t_star(7600);
t_star(1:7600) = 0;

%% Ramp kinematics
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

%  writerObj = VideoWriter('fml_3000_EoS.avi');
%  open(writerObj)

%% Generate Pressure Map

x = [0.1 0.1 0.1 0.1 0.2 0.2 0.5 0.2 0.5 0.3 0.3 0.4 0.3 0.4];  
y = [0.2 0.4 0.6 0.8 0.4 0.6 0.6 0.8 0.8 0.4 0.6 0.6 0.8 0.8];

for i = 8100 %7701:8:9000 %5501:8:10080%:10:10300 %8350 9850 9200 8281
    
tStep1 = 5500; %find(s_star > 10,1);
tStep2 = 7600;%find(s_star > 18,1);

% p1 = mean(Cp1(tStep1:tStep2,:));
% p2 = mean(Cp2(tStep1:tStep2,:));
% p3 = mean(Cp3(tStep1:tStep2,:));
% p4 = mean(Cp4(tStep1:tStep2,:));

if i > 7840
d = SStr{1,1}(1,i-7840) - SStr{1,1}(1,i-7839);
p1 = Cp1(i,:)./((d*c/T_s)^2);
p2 = Cp2(i,:)./((d*c/T_s)^2);
p3 = Cp3(i,:)./((d*c/T_s)^2);
p4 = Cp4(i,:)./((d*c/T_s)^2);
else
p1 = Cp1(i,:);
p2 = Cp2(i,:);
p3 = Cp3(i,:);
p4 = Cp4(i,:);
end

p1 = [p1(:,1:4),p1(:,6:end)];
p2 = [p2(:,1:4),p2(:,6:end)];
p3 = [p3(:,1:4),p3(:,6:end)];
p4 = [p4(:,1:4),p4(:,6:end)];

% for i = 1:length(Cp1(1,:))
% p1(i) = pSample(i);
% end
             
[xG, yG] = meshgrid(-0.01:0.001:1, -0.01:0.001:1);

vG(:,:,1) = griddata(x,y,p1,xG,yG,'v4');
vG(:,:,2) = griddata(x,y,p2,xG,yG,'v4');
vG(:,:,3) = griddata(x,y,p3,xG,yG,'v4');
vG(:,:,4) = griddata(x,y,p4,xG,yG,'v4');

%% Colormap

l=50;
n=50;
m=10;

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
%plotPressureMaps(xG, yG, vG, x, y, P, [-1.5 0],i,t_star)

%% Plot Map
 fig=figure(1);
    set(fig,'position',[0.5 0.5 11 7],'color',[1 1 1],'visible','on')
    set(gca,'fontsize',18)

    xx = [0  1 -1];
    yy = [0  1  1];
    Q = patch(xx,yy,[0.75 0.75 0.75]);
    Q.FaceColor = [0.75 0.75 0.75];
    Q.EdgeColor = [0 0 0];
    Q.LineWidth = 2;

    hold on
    % h = contourf(Xgr,Ygr,Vgr,linspace(-1,1,200),'LineStyle','none');
    h = surf(xG,yG,vG(:,:,3));
    h.EdgeColor = 'none';
    caxis([-2 0])

    % cmap = four_color_vort;
    % colormap(cmap(1:112,:)); %only uses the cold portion of the color map
    %colormap(cmap(112:225,:)); %only uses the warmportion of the color map
    colormap('jet')

    %colorbar
    hold on
    l = plot3(x,y,P(3,:),'.k','markersize',10);
    axis equal
    view(90,-90)
    %colorbar
     xlabel('$z/c$','interpreter','Latex','fontsize',36)
     ylabel('$x/c$','interpreter','Latex','fontsize',36)
    %title('$\alpha = 20^{\circ}$, $s^*=1$','interpreter','Latex','fontsize',30)
    set(gca,'FontSize',22)
    axis([-0.1 1 -0.2 1 -10 10])
    set(gca,'XTick',[-0.1 0:0.2:1])
    set(gca,'YTick',[-0.2:0.2:1])
    ax = gca;
    ax.GridLineStyle = '--';
    ax.FontSize = 28;
    ax.GridAlpha = 0.5;
    ax.LineWidth = 2;
    grid on
    
    if i < 7840
    t = annotation('textbox');
    t.String = 'Steady State';
    t.Interpreter = 'latex';
    t.FontSize = 30;
    t.LineStyle = 'none';
    t.Position = [0.3 0.8 0.1 0.1];
    else
    str = num2str(SStr{1,1}(1,i-7840),'%0.2f');
    t = annotation('textbox');
    %t.String = ['$s^*$ =', str];
    t.String = ['$S^* =$',str];
    t.Interpreter = 'latex';
    t.FontSize = 30;
    t.LineStyle = 'none';
    t.Position = [0.3 0.8 0.1 0.1];
    end
     grid off
     
    B=colorbar; 
    set(B, 'Position', [0.83 .11 .06 .85])
    set(get(B,'ylabel'),'String','$C_p$','Interpreter','Latex','FontSize', 30);
     cpos = B.Position;
     cpos(3) = 0.5*cpos(3);
     B.Position = cpos;
 

    %% Make Movie
%      N = myaa;
%     temp =getframe(gcf);
%     writeVideo(writerObj,temp);
%     disp(i)
%     close all
    
end

%     close(writerObj);

%% Plot pressure traces
% t_star = t*U/c;
% 
% fig=figure(2);
% set(fig,'position',[1 1 8 6],'color',[1 1 1])
% set(gca,'fontsize',13)
% 
% plot(s{1,1}(1,:), Cp1(1162:L,3), 'c', ...
%      s{2,1}(1,:), Cp2(1122:L,3), 'm', ...
%      s{1,1}(1,:), Cp3(1162:L,3), 'g', ...
%      s{2,1}(1,:), Cp4(1122:L,3), 'r',  'linewidth',1.35);
%  
% xlabel('$s^*$','Interpreter','Latex','FontSize', 17);
% ylabel('$C_p$','Interpreter','Latex','FontSize',17);
% axis([0 40 -4 1])
% grid on
% ax = gca;
% ax.GridLineStyle = '--';
% 
% ax = gca;
% ax.GridLineStyle = '--';
% ax.LineWidth = 2;
% set(gca,'position',[0.65 0.6 6.4 4.6])
% legend('20degStar1','20degStar2','30degStar1','30degStar2')
% hold on
% grid on
   
%%

% s_star = s_star';
% 
% subplot(2,2,1)
% plot(t_star, Cp1(:,3), 'c', ...
%      t_star, Cp2(:,3), 'm', ...
%      t_star, Cp3(:,3), 'g', ...
%      t_star, Cp4(:,3), 'r',  'linewidth',1.35);
%  
% xlabel('$s^*$','Interpreter','Latex','FontSize', 17);
% ylabel('$C_p$','Interpreter','Latex','FontSize',17);
% axis([0 40 -4 1])
% grid on
% ax = gca;
% ax.GridLineStyle = '--';

% subplot(2,2,2)
% plot(t_star, Cp1(:,7), 'c', ...
%      t_star, Cp2(:,7), 'm', ...
%      t_star, Cp3(:,7), 'g', ...
%      t_star, Cp4(:,7), 'r',  'linewidth',1.35);
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

%% Trill shit dawg.. seriously don't bother with this crap here

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
%     t4 = dSSR/Uf;
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
% 
