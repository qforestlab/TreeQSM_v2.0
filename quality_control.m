%QSM quality control
%Branch by Zane 17/05/23, made for overlay of point cloud and cylinders
%Last update 08/05/24 to make cylinder image display more similar to TreeQSM v2.4
clc;
clear all;
cd
clouds = dir('pointcloud\*.txt');
models_dir = ['results'];
for i = 1:length(clouds)
    P = dlmread([clouds(i).folder,'\',clouds(i).name]);
    pc = P(:,1:3);

    %normalize x/y/z values for graphing
    %ncols = size(pc)
    %ncols = ncols(2)
    %for j = 1:ncols
    %    pc(:,j) = pc(:,j) - min(pc(:,j))
    %end

    C = strsplit(clouds(i).name,{'.txt'}); 
    name = dir(['results/',C{1},'-*.mat']);
    QSMs(1) = load(name.name);
    model = QSMs(1);

    disp(name)
    plot_cylinder_model(model.Rad,model.Len,model.Axe,model.Sta,model.BoC,1,4,0.75);
    figname_cyl = ['figures\',C{1},'_cylinders_v2.0.png'];
    az=0;
    el=0;
    grid off
    view (az,el)
    set(gca, 'XTickLabel',{},'ZTickLabel',{});
    axis off
    set(gcf, 'Position');
    saveas(gcf,figname_cyl)
    %option for high resolution export
    exportgraphics(gcf,figname_cyl,'Resolution',1000)

    hold on
    plot_point_cloud(pc,1,1,'k');
    hold off
    set(gcf, 'Position', get(0, 'Screensize'));
    az=0;
    el=0;
    grid off
    view (az,el)
    figname_combine = ['figures\',C{1},'_combine_v2.0.png'];
    saveas(gcf,figname_combine)
    pause;
end