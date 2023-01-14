%% Plotting pressure distribution over NACA0012 delta wing
%  and make animation videos
%
%
%     Required functions:  PlotWing.m; linearKinematics.m; smoother200.m
%     Required files:      colorbar.mat; Length.mat

close all; clear all; clc;      %#ok<*CLALL>

close all; clear all; clc;
%set(0,'DefaultFigureUnits','Inches');
%set(0,'DefaultAxesUnits','Inches');
%set(0, 'DefaultFigurePosition', [14.4375   -2.3958  10.0000    6.0000])

%% Load Data

load('dataset_allCases_inOrder_P.mat')
L = 2500; % length of each case for ML, 2500 by default, including 2000 points
          % after gust and 500 points before gust
          
caseToPlot = 16; % which case to plot, just for single contour, not animation

%% Input Variables
U = 1;
Uf = 1.5;
c = 0.3;
T = 0.3; % s periodic of the gust = 0.3s
sampleRate = 1000;
T_s = 1/sampleRate;
t = linspace(0,T_s*L,L);
t_star = t/T;

%% Movie

writerObj = VideoWriter(['contour_case' num2str(caseToPlot,'%02i') '.avi']);
open(writerObj);

%% Generate Pressure Map

x = [ 0.1 0.1 0.2 0.3 0.1 0.2 0.3 0.4 0.5 0.1 0.2 0.3 0.4 0.5 0.7];
y = [ 0.2 0.4 0.4 0.4 0.6 0.6 0.6 0.6 0.6 0.8 0.8 0.8 0.8 0.8 0.8];

%%
for i = 300:5:1300 % 300:5:1300 for movie % 680(case21) 890(case17) (if single number, just for figure, not for movie)

    data_P_temp = data_P((L*(caseToPlot-1)+1):(L*caseToPlot),1:15);

    p1 = data_P_temp(i,1:15); % port 16 links to stagnation point, so no use port 16

    chord = 1;
    grSize = 500;
    [xx,yy]=meshgrid(linspace(-1,1,grSize),linspace(0,1,grSize));

    vG(:,:,1) = griddata(x,y,p1,xx,yy,'natural');


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

    P = [p1; p1; p1; p1];

    plotWing(xx,yy,vG,grSize,1,[-3 3],i)  % [-7 2] for all cases, for consistant
    % [-2 1] case01; [-3 1.5] case02; [-5 1] case03; [-7 2] case04; [-2 1] case05 
    % [-3 2] case06; [-3.5 0] case07; [-2 0] case08; [-4 2] case09; 
    % [ -3.5 1] case10; [-2 2] case11; [-3 4] case12; [-8 2] case13;
    % [-6 1] case14; [-3.5 1] case15; [-3 3] case16
    % [-3 1] case21
    % plotPressureMaps(xx, yy, vG, x, y, P, [-1.5 0],i,t_star,SStr)
    magma=magma();

        R = [zeros(1,300),linspace(0,1),ones(1,100),linspace(1,0.419608), ...
        linspace(0.419608,1)];
        G = [zeros(1,100),linspace(0,1),linspace(1,0.501961), ...
            linspace(0.501961,1),linspace(1,0.380392),linspace(0.380392,0) ... 
            linspace(0,1)];
        B = [linspace(0.392157,1),ones(1,100),linspace(1,0),zeros(1,300), ...
            linspace(0,1)];

        cpDesat = [R',G',B'];
        colormap(cpDesat)

        if i < 500
            str = num2str(-(500-i)*T_s/T,'%0.2f');
            text = annotation('textbox');
            text.String = ['$T^*$ = ', str];
            text.Interpreter = 'latex';
            text.FontSize = 32;
            text.LineStyle = 'none';
            text.Position = [0.05 0.75 0.1 0.1];
            Color20 = [1, 0.25, 0];
            Color30 = [0,0.0,1];
        else
            str = num2str((i-500)*T_s/T,'%0.2f');
            text = annotation('textbox');
            text.String = ['$T^*$ = ', str];
            text.Interpreter = 'latex';
            text.FontSize = 32;
            text.LineStyle = 'none';
            text.Position = [0.05 0.75 0.1 0.1];
            Color20 = [1, 0.25, 0];
            Color30 = [0,0.0,1];
        end
        
        % Make Movie
        temp =getframe(gcf);
        writeVideo(writerObj,temp);
        disp(i)
        close all
        
end
% %%
close(writerObj);

