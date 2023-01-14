function [t,b] = plotSlicedWingY(theta,xCut)


chord = 1;
grSize = 1000;
[xx,yy]=meshgrid(linspace(-chord,chord,grSize),linspace(0,chord,grSize));
c=[linspace(0,1,grSize/2),linspace(1,0,grSize/2)]*chord;

zz = zeros(size(xx));
k = 1;

for ii=1:grSize
    zz(1:sum(yy(:,ii)<=c(ii)),ii)=flipud(5*0.12*c(ii)*(...
        0.2969*(yy(1:sum(yy(:,ii)<=c(ii)),ii)./c(ii)).^(1/2)...
        -0.1260*(yy(1:sum(yy(:,ii)<=c(ii)),ii)./c(ii)).^(1)...
        -0.3516*(yy(1:sum(yy(:,ii)<=c(ii)),ii)./c(ii)).^(2)...
        +0.2843*(yy(1:sum(yy(:,ii)<=c(ii)),ii)./c(ii)).^(3)...
        -0.1015*(yy(1:sum(yy(:,ii)<=c(ii)),ii)./c(ii)).^(4)));    
    
end

%remove excess
zz(1,1) = 0;

for i = 1:grSize
  j = 1;
  
  while zz(j,i) ~= 0 
      j = j+1;
  end
  
  if i == grSize/2 + 2 || i == grSize/2 + 1
    xx(j+4:end,i)= nan;
    yy(j+4:end,i)= nan;
    zz(j+4:end,i)= nan;
  elseif i == grSize/2 - 2 || i == grSize/2 - 1
    xx(j+4:end,i)= nan;
    yy(j+4:end,i)= nan;
    zz(j+4:end,i)= nan;
  else
    xx(j+3:end,i)= nan;
    yy(j+3:end,i)= nan;
    zz(j+3:end,i)= nan;
  end
  
end



lx = linspace(-0.3,0,grSize);
ly = linspace(0,0.3,grSize);
lxN = linspace(0,0.3,grSize);
lyN = linspace(0.3,0,grSize);

% for j = 1:length(xx);
%     xx(j,round(j/2)) =  lx(j);
%     yy(j,round(j/2)) =  ly(j);
%     zz(j,round(j/2)) = 0;
%     xx(j,grSize-round(j/2)) = lxN(j);
% end


zzB = -zz;
zz = zz - 0.07*chord;
zzB = zzB  - 0.07*chord;

%% Rotate translate

yOffset = chord/2; 
yy = yy - yOffset;

yy = yy.*cosd(theta) - zz.*sind(theta);
zz = yy.*sind(theta) + zz.*cosd(theta);
zzB = yy.*sind(theta) + zzB.*cosd(theta);

yy = yy + yOffset;
yy = yy + 1.2*chord;
xx = xx + 0.02;
Cut = xCut/chord;



if xCut > 0
jj = Cut*grSize;
yy(jj:end,:) = nan;


t =  surf(xx,zz,yy);
t.FaceColor = [0.75,0.75,0.75];
t.LineStyle = 'none';
hold on
b =  surf(xx,zzB,yy);
b.FaceColor = [0.75,0.75,0.75];
b.LineStyle = 'none';
else
    t = 0;
    b = 0;
end

end

