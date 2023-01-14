function plotPressureMaps(Xgr,Ygr,Vgr,XX,YY,Pres,Range,idx,t_star,s)

    fig=figure(1);
    set(fig,'position',[0.1 0.1 10 6],'color',[1 1 1],'visible','on')
    set(gca,'fontsize',18)

    %% Plot 1
    subplot(2,2,1)

    xx = [0  1 -1];
    yy = [0  1  1];
    Q = patch(xx,yy,[0.75 0.75 0.75]);
    Q.FaceColor = [0.75 0.75 0.75];
    Q.EdgeColor = 'none';

    hold on
    % h = contourf(Xgr,Ygr,Vgr,linspace(-1,1,200),'LineStyle','none');
    h = surf(Xgr,Ygr,Vgr(:,:,1));
    h.EdgeColor = 'none';
    caxis(Range)

    % cmap = four_color_vort;
    % colormap(cmap(1:112,:)); %only uses the cold portion of the color map
    %colormap(cmap(112:225,:)); %only uses the warmportion of the color map
    colormap('jet')

    %colorbar
    hold on
    l = plot3(XX,YY,Pres(1,:),'.k','markersize',10);
    axis equal
    view(90,-90)
    %colorbar
    xlabel('$z/c$','interpreter','Latex','fontsize',16)
    %ylabel('$x/c$','interpreter','Latex','fontsize',16)
    title('$\alpha = 20^{\circ}$, $g^*=1$','interpreter','Latex','fontsize',16)
    set(gca,'FontSize',16)
    axis([-0.2 1 -0.2 1 -10 10])
%     colormap off
    colormap default
    
%     if idx < 7600
%     t1 = annotation('textbox');
%     t1.String = 'Steady State';
%     t1.Interpreter = 'latex';
%     t1.FontSize = 23;
%     t1.LineStyle = 'none';
%     t1.Position = [0.6 0.2 0.1 0.1];
%     else
%     str = num2str(t_star(idx),'%0.2f');
%     t1 = annotation('textbox');
%     t1.String = ['$t^*$ =', str];
%     t1.Interpreter = 'latex';
%     t1.FontSize = 23;
%     t1.LineStyle = 'none';
%     t1.Position = [0.6 0.2 0.1 0.1];
%     end
%     
    %% Plot 2
    subplot(2,2,2)

    xx = [0  1 -1];
    yy = [0  1  1];
    Q = patch(xx,yy,[0.75 0.75 0.75]);
    Q.FaceColor = [0.75 0.75 0.75];
    Q.EdgeColor = 'none';

    hold on
    % h = contourf(Xgr,Ygr,Vgr,linspace(-1,1,200),'LineStyle','none');
    h = surf(Xgr,Ygr,Vgr(:,:,3));
    h.EdgeColor = 'none';
    caxis(Range)

    % cmap = four_color_vort;
    % colormap(cmap(1:112,:)); %only uses the cold portion of the color map
    %colormap(cmap(112:225,:)); %only uses the warmportion of the color map
    colormap('jet')

    %colorbar
    hold on
    l = plot3(XX,YY,Pres(3,:),'.k','markersize',10);
    axis equal
    view(90,-90)
    %colorbar
    xlabel('$z/c$','interpreter','Latex','fontsize',16)
    %ylabel('$x/c$','interpreter','Latex','fontsize',16)
    title('$\alpha = 30^{\circ}$, $g^*=1$','interpreter','Latex','fontsize',16)
    set(gca,'FontSize',16)
      axis([-0.2 1 -0.2 1 -10 10])
      
    %% Plot 4
    subplot(2,2,3)

    xx = [0  1 -1];
    yy = [0  1  1];
    Q = patch(xx,yy,[0.75 0.75 0.75]);
    Q.FaceColor = [0.75 0.75 0.75];
    Q.EdgeColor = 'none';

    hold on
    % h = contourf(Xgr,Ygr,Vgr,linspace(-1,1,200),'LineStyle','none');
    h = surf(Xgr,Ygr,Vgr(:,:,2));
    h.EdgeColor = 'none';
    caxis(Range)

    % cmap = four_color_vort;
    % colormap(cmap(1:112,:)); %only uses the cold portion of the color map
    %colormap(cmap(112:225,:)); %only uses the warmportion of the color map
    colormap('jet')

    %colorbar
    hold on
    l = plot3(XX,YY,Pres(2,:),'.k','markersize',10);
    axis equal
    view(90,-90)
    %colorbar
    xlabel('$z/c$','interpreter','Latex','fontsize',16)
    ylabel('$x/c$','interpreter','Latex','fontsize',16)
    title('$\alpha = 20^{\circ}$, $g^*=2$','interpreter','Latex','fontsize',16)
    set(gca,'FontSize',16)
     axis([-0.2 1 -0.2 1 -10 10])
    %% Plot 4
      subplot(2,2,4)

    xx = [0  1 -1];
    yy = [0  1  1];
    Q = patch(xx,yy,[0.75 0.75 0.75]);
    Q.FaceColor = [0.75 0.75 0.75];
    Q.EdgeColor = 'none';

    hold on
    % h = contourf(Xgr,Ygr,Vgr,linspace(-1,1,200),'LineStyle','none');
    h = surf(Xgr,Ygr,Vgr(:,:,4));
    h.EdgeColor = 'none';
    caxis(Range)

    % cmap = four_color_vort;
%      colormap(cmap(1:112,:)); %only uses the cold portion of the color map
    %colormap(cmap(112:225,:)); %only uses the warmportion of the color map
    colormap('jet')

    %colorbar
    hold on
    l = plot3(XX,YY,Pres(4,:),'.k','markersize',10);
    axis equal
    view(90,-90)
    %colorbar
    xlabel('$z/c$','interpreter','Latex','fontsize',16)
    title('$\alpha = 30^{\circ}$, $g^*=2$','interpreter','Latex','fontsize',16)
    ylabel('$x/c$','interpreter','Latex','fontsize',16)
    axis([-0.2 1 -0.2 1 -10 10])
     set(gca,'FontSize',16)
     
    if idx < 7840
    t = annotation('textbox');
    t.String = 'Steady State';
    t.Interpreter = 'latex';
    t.FontSize = 16;
    t.LineStyle = 'none';
    t.Position = [0.45 0.5 0.1 0.1];
    else
    str = num2str(s{1,1}(1,idx-7840),'%0.2f');
    t = annotation('textbox');
    t.String = ['$s^*$ =', str];
    t.Interpreter = 'latex';
    t.FontSize = 16;
    t.LineStyle = 'none';
    t.Position = [0.45 0.5 0.1 0.1];
    end
    
    B=colorbar; 
    set(B, 'Position', [0.9 .11 .04 .8150])
    set(get(B,'ylabel'),'String','$C_p$','Interpreter','Latex','FontSize', 20);
    cpos = B.Position;
    cpos(3) = 0.5*cpos(3);
    B.Position = cpos;
    
%colorbar off
end