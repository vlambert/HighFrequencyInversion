function  MakeFrame(tmp,tmp2,x_bp,y_bp,dx,dy,t,x_ev,y_ev,frameDir,ii,nt)
    h1 = figure('Visible','Off');clf
    set(gcf,'Position',[1 1 996 384]);
    set(gcf,'color','w');
    subplot(3,2,[1 3 5]);
    h=pcolor(x_bp-dx/2,y_bp-dy/2,tmp');
    set(h,'EdgeColor','none');
    axis equal;
    xlim([min(x_bp)-dx/2 max(x_bp)-dx/2])
    ylim([min(y_bp)-dy/2 max(y_bp-dy/2)])
    caxis([0 1]);
    if ii == (nt+1)
        title('Linear Stack');
    else
        title(sprintf('Linear Stack: Time %.2f',t(ii)));
    end
    hold on; plot(x_ev,y_ev,'rp','MarkerSize',10);
    %plot(xpeak(1:ii),ypeak(1:ii),'bs');
    
    % Track the temporal progress along the sample grid point
    %subplot(4,2,7);hold on;
    %if(ii<nt+1) 
    %    plot(t(ii),sample1(ii),'rx');
    %end

    
    %subplot(4,2,[2 4 6]);hold off;
    subplot(3,2,[2 4 6]);
    h=pcolor(x_bp-dx/2,y_bp-dy/2,tmp2');
    set(h,'EdgeColor','none');
    axis equal;
    xlim([min(x_bp)-dx/2 max(x_bp)-dx/2])
    ylim([min(y_bp)-dy/2 max(y_bp)-dy/2])
    caxis([0 1]);
    title('Square Stack');
    hold on; plot(x_ev,y_ev,'rp','MarkerSize',10);
    %plot(xpeaks(1:ii),ypeaks(1:ii),'bs');
    %subplot(4,2,8);hold on;
    %if(ii<nt+1)
    %    plot(t(ii),sample2(ii),'rx');
    %end
    %mov(ii)=getframe(gcf);
    %img2=getframe(gcf);
    %imwrite(img2.cdata, [frameDir,sprintf('Frames_%d.png', ii)]);
    saveas(gcf,[frameDir,sprintf('Frames_%d', ii)],'png')

end