%Andrew Burt - a.burt.12@ucl.ac.uk

function plotCyls(infile)
    h = figure('visible','off');
    %campos([-20,18,8]);
    view([0 30])
    hold on;
    axis equal;
    grid off;
    axis off;
    for z=1:length(infile)
        c = rand(1,3);
        load(infile{z});
        for i = 1:length(Rad)
            [X,Y,Z] = cylinder2P(Rad(i),8,[Sta(i,1),Sta(i,2),Sta(i,3)],[Sta(i,1)+(Len(i)*Axe(i,1)),Sta(i,2)+(Len(i)*Axe(i,2)),Sta(i,3)+(Len(i)*Axe(i,3))],c);
        end
    end
    ti = get(gca,'TightInset');
    set(gca,'Position',[ti(1) ti(2) 1-ti(3)-ti(1) 1-ti(4)-ti(2)]);
    set(gca,'units','centimeters')
    pos = get(gca,'Position');
    ti = get(gca,'TightInset');
    set(gcf, 'PaperUnits','centimeters');
    set(gcf, 'PaperSize', [pos(3)+ti(1)+ti(3) pos(4)+ti(2)+ti(4)]);
    set(gcf, 'PaperPositionMode', 'manual');
    set(gcf, 'PaperPosition',[0 0 pos(3)+ti(1)+ti(3) pos(4)+ti(2)+ti(4)]);
    print(h,'-dpdf','-r600','models.pdf');
end
