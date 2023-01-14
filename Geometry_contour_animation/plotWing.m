function plotWing(XX,YY,VG,gridSz,pSet,range,idx)

chord = 1;
grSize = gridSz;
if pSet == 1 || pSet == 2
    theta = 20;
elseif pSet == 3 || pSet == 4
    theta = 30;
else
    fprintf('youGoofed')
end

c=[linspace(0,1,grSize/2),linspace(1,0,grSize/2)]*chord;
zz = zeros(size(XX));

for ii=1:grSize
    zz(1:sum(YY(:,ii)<=c(ii)),ii)=flipud(5*0.12*c(ii)*(...
        0.2969*(YY(1:sum(YY(:,ii)<=c(ii)),ii)./c(ii)).^(1/2)...
        -0.1260*(YY(1:sum(YY(:,ii)<=c(ii)),ii)./c(ii)).^(1)...
        -0.3516*(YY(1:sum(YY(:,ii)<=c(ii)),ii)./c(ii)).^(2)...
        +0.2843*(YY(1:sum(YY(:,ii)<=c(ii)),ii)./c(ii)).^(3)...
        -0.1015*(YY(1:sum(YY(:,ii)<=c(ii)),ii)./c(ii)).^(4)));    
    
end

%remove excess
zz(1,1) = 0;

for i = 1:grSize
  j = 1;
  
  while zz(j,i) ~= 0 
      j = j+1;
  end
  
  if i == grSize/2 + 2 || i == grSize/2 + 1
    XX(j+4:end,i)= nan;
    YY(j+4:end,i)= nan;
    zz(j+4:end,i)= nan;
  elseif i == grSize/2 - 2 || i == grSize/2 - 1
    XX(j+4:end,i)= nan;
    YY(j+4:end,i)= nan;
    zz(j+4:end,i)= nan;
  else
    XX(j+3:end,i)= nan;
    YY(j+3:end,i)= nan;
    zz(j+3:end,i)= nan;
  end
  
end

lx = linspace(-0.3,0,grSize);
ly = linspace(0,0.3,grSize);
lxN = linspace(0,0.3,grSize);
lyN = linspace(0.3,0,grSize);

zzB = -zz;
zz = zz - 0.07*chord;
zzB = zzB  - 0.07*chord;

%% Rotate translate

yOffset = chord/2; 
YY = YY - yOffset;

YY = YY.*cosd(theta) - zz.*sind(theta);
zz = YY.*sind(theta) + zz.*cosd(theta);
zzB = YY.*sind(theta) + zzB.*cosd(theta);

YY = YY + yOffset;
YY = YY + 1.2*chord;
XX = XX + 0.02;


fig=figure(1);
    set(fig,'position',[200 200 1000 500],'color',[1 1 1],'visible','on')
    set(gca,'fontsize',18)
    
t =  surf(XX,YY,zz);
t.FaceColor = [0.75,0.75,0.75];
t.LineStyle = 'none';
hold on
b =  surf(XX,YY,zzB);
b.FaceColor = [0.75,0.75,0.75];
b.LineStyle = 'none';

tmp = flipud(VG(:,:,pSet));
VG = fliplr(tmp);
f = surf(XX,YY,zz+0.01,VG);
f.LineStyle = 'none';
caxis(range)
colormap('jet')
ax = gca;
%set(ax, 'Visible','off')
axis off

axis equal
xlim([-1 1])
cameratoolbar
cameratoolbar('ResetCamera')
camproj('perspective')
cameratoolbar('SetCoordSys','z')
% camva(3.4120)
% camup([0.0517    0.1709    0.9839])
% campos([-2.9707   -8.6486    1.9066])
% camtarget([0.1486    1.6628   -0.0483])
% camva(3.4120)
% camup([ 0.1008    0.2086    0.9728])
% campos([-4.8938  -16.0721    3.3646])
% camtarget([0.4475    1.5848    0.0171])
% camva(3.4120)
% camup([0.1008    0.2086    0.9728])
% campos([-7.5025  -14.8397    4.3255])
% camtarget([  0.4364    1.5802   -0.0177])
% light('Position',[2 2  6],'Style','infinite')
% lighting gouraud
% material dull
% camva(3.4120)
% camup([0.1008    0.2086    0.9728])
% campos([-5.2786   -9.6841    2.9586])
% camtarget([0.2155    1.6792   -0.0471])
camva(3.7385)
camup([0.1008    0.2086    0.9728])
campos([-5.2786   -9.6841    2.9586])
camtarget([0.2155    1.6792   -0.0471])
light('Position',[2 2  6],'Style','infinite')
lighting gouraud
material dull


   % B=colorbar;
%     set(B, 'Position', [0.11 0.1 0.04 0.8])
    % set(get(B,'ylabel'),'String','$C_p$','Interpreter','Latex','FontSize', 32);
%     cpos = B.Position;
%     cpos(3) = 0.5*cpos(3);
%     B.Position = cpos;
B=colorbar;
%set(get(B,'ylabel'),'String','$\omega_z^*$','Interpreter','Latex','FontSize', 45);
set(get(B,'ylabel'),'String','$C_p$','Interpreter','Latex','FontSize', 45);
set(B, 'Position', [0.85 0.125 0.04 0.8])
    cpos = B.Position;
    cpos(3) = 0.5*cpos(3);
    B.Position = cpos;
   set(gca,'fontsize', 32);
   set(gca,'fontsize', 30);
   %colorbar 
   
fprintf('plotted')

end


