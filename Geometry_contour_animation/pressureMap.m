close all; clear all; clc;
set(0,'DefaultFigureUnits','Inches');
set(0,'DefaultAxesUnits','Inches');

%% Input Variables
sampleRate = 1000;
T_s=1/sampleRate;
t =    linspace(0,T_s*L,L);
t_star = t*U/c;
s_star = t_star/2;


%% Generate Pressure Map

x = [0.1 0.3 0.5 0.7 0.1 0.2 0.3 0.4 0.1 0.1
    
        ];
    
y = [0.2 0.4 0.6 0.8 0.6 0.6 0.6 0.6 0.4 0.8];

load('meanLE.mat')
load('meanSpanwise.mat')
load('meanChord.mat')

% P10 = -1*[1.229 1.083 (1.133+1.018)/2 0.9983 0.4262 0.4823 0.5845 0.6788 ...
%           0.4797 0.1435];
% 
% P20 = -1*[1.4479 1.1684 0.8783 0.5871 0.5921 0.6451 0.9084 0.9064  ...
%           0.9587 0];
% 
% P30 = -1*[1.1348 0.9346 0.7721 0.5533 0.8051 0.8716 0.9235 0.8762  ...
%            1.161 1.347];

for i = 1:3
p{i} = [(leadingEdge(i,1)+chord(i,1))/2 leadingEdge(i,2) ...
        (leadingEdge(i,3)+spanWise(i,5))/2 leadingEdge(i,4) ...
         (spanWise(i,1)+chord(i,3))/2 spanWise(i,2) spanWise(i,3) ...
         spanWise(i,4) chord(i,2) chord(i,4)];
end
             
[xG, yG] = meshgrid(-0.01:0.001:1, -0.01:0.001:1);

vG = griddata(x,y,p{3},xG,yG,'cubic');

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

%% Plot pressure map
fig=figure(1);
set(fig,'position',[1.5 1.5 10 6],'color',[1 1 1])
set(gca,'fontsize',18)

xx = [0  1 -1];
yy = [0  1  1];
P = patch(xx,yy,[0.75 0.75 0.75]);
P.FaceColor = [0.75 0.75 0.75];
P.EdgeColor = 'none';

hold on
% h = contourf(xG,yG,vG,linspace(-1,1,200),'LineStyle','none');
h = surf(xG,yG,vG);
h.EdgeColor = 'none';
caxis([-1 0])
    
load('colorbar.mat')
cmap = four_color_vort;
colormap(cmap(1:112,:)); %only uses the cold portion of the color map
%colormap(cmap(112:225,:)); %only uses the warmportion of the color map
%colormap('jet')

%colorbar
hold on
l = plot3(x,y,p{3},'.k','markersize',10);
axis equal
view(180,-90)

xlabel('$z/c$','interpreter','Latex','fontsize',30)
ylabel('$x/c$','interpreter','Latex','fontsize',30)
