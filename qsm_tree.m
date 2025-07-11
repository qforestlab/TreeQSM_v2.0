function [TreeData, Sta, Axe, Rad, Len, BoC, CPar, CExt, BOrd, BPar, ...
    BVol, BLen, BAng, BSeg, FCB, BChi, CiB, CChi, CiS, Added, P, Bal, ...
    Cen, Nei, Segs, SPar, SChi, SoC, Base, Forb]...
    = qsm_tree(P,dmin0,rcov0,nmin0,dmin,rcov,nmin,lcyl,NoGround,string, ...
        rfil1, nfil1, rfil2, nfil2)

% ---------------------------------------------------------------------
% QSM_TREE.M     Creates quantitative structure tree models from point 
%                   clouds scanned from trees.
%
% Version 2.0
% Author        Pasi Raumonen
% Created       15 March 2014
%
% This software cannot be distributed without the permission of P. Raumonen
% and it is only for non-commercial use. These restrictions holds also for 
% the modules or subprograms of the software:
% FILTERING.M, COVER_SETS.M, TREE_SETS.M, SEGMENTS.M, CORRECT_SEGMENTS.M,
% RELATIVE_SIZE.M, CYLINDERS.M, BRANCHES.M, TREE_DATA.M
% ---------------------------------------------------------------------

% Produces a cylindrical tree model from point cloud, which is a sample of
% the tree surface. Uses covers of point cloud for segmentation of the
% point cloud into stem and branches. The cover sets are with small sets, 
% which are along the surface, such that each point belongs at most one 
% cover set; i.e. the cover is a partition of the point cloud. The
% segmentation process uses two different covers: First cover is with large
% size sets and the purpose of the first cover is to separate the
% measurements belonging to the tree from those of ground and understorey
% measurements. The second purpose is to make a rough segmentation of the
% tree sets into stem and branches so that we can generate a new finer
% cover that also refines further along the branches, branch order and 
% height. That is the first cover may have size e.g. dmin0 = 8cm where as 
% the second cover has the size at the base of the stem e.g. dmin = 2cm and 
% then half the size, i.e. 1cm, at the top of the tree and at the tip of the
% each branch. Furthermore, the segmentation is corrected so that first the
% stem is modified to be from the base to top of the tree. Then each branch
% (in the increasing branching order) is modified similarly so that it is
% from its base to the furthest tip of its child segments. When the second
% finer cover is generated and segmented and corrected into stem and
% branches, then the cylinder model is generated by fitting cylinders into
% subsegments or regions of each branch/segment. The length of these
% regions is controlled so that the relative length of the cylinder, lcyl =
% length/radius is about equal to the given value. During the cylinder
% fitting process the radius of the cylinders is controlled so that it is
% not too large or small and this uses the size of the parent cylinder and
% parabola fitted into the branch taper.
% 
% The cover is generated such that at first the point cloud is covered 
% with balls with radius RCOV. This first cover is such that 
% 1) the minimum distance between the centers is DMIN, and 
% 2) the maximum distance from any point to nearest center is also DMIN.
% Then the first cover of RCOV-balls is used to define a second cover:
% each RCOV-ball A defines corresponding cover set B in the second cover
% such that B contains those points of A that are nearer to the center of
% A than any other center of RCOV-balls. The RCOV-balls also define 
% the neighbors for the second cover: Let CA and CB denote cover sets in 
% the and second cover, and BA and BB their RCOV-balls. Then CB is 
% a neighbor of CA, and vice versa, if BA and CB intersect or 
% BC and CA intersect. There is, however, exception to this rule: if the
% gap between the sets is too large, then they are not neighbors.
%
% The point cloud is first filtered, if RFIL1 > 0, to remove isolated 
% points and point groups. This is done as follows: First RFIL1-ball is
% generated for each point and then if a ball contains less than NFIL1 
% points, the point is removed. Next the remaining point cloud is covered 
% with RFIL2-balls, and the connected components are determined. If 
% a component has less than NFIL2 cover sets, then the points in the 
% component is removed.
%
% Inputs: 
% P         (Un)filtered point cloud, (m_points x 3)-matrix, the rows
%               give the coordinates of the points.
%               The order of the points is not meaningful
% dmin      Minimum distance between centers of cover sets; i.e. the
%               minimum diameter of a cover set
% rcov      Radius of the balls used to generate the cover sets, these 
%               balls are also used to determine the neighbors and the 
%               cover set characteristics              
% nmin      Minimum number of points in a rcov-ball
% lcyl      Cylinder length/radius ratio
% NoGround  Logical value, true if no ground in the point cloud, in which
%               case defines the base of the trunk as the lowest part the
%               cloud
% string    Name string for saving output files
% rfil1     Radius of cover sets used in the first filtering process
% nfil1     Minimum number of points in the cover sets passing the first filtering
% rfil2     Radius of cover sets used in the second filtering process
% nfil2     Minimum number of cover sets in components passing the second filtering
% 
% Cylinder model outputs:
% Sta       Starting points of the cylinders, matrix
% Axe       Axes of the cylinders, matrix
% Rad       Radii of the cylinders, vector
% Len       Lengths of the cylinders, vector
% CPar      Parent cylinder of each cylinder, vector
% CExt      Extension cylinder of each cylinder, vector
% BoC       Branch of the cylinder, vector
% BOrd      Branch order, vector
% BPar      Parent branch, vector
% BVol      Volumes of the branches, vector
% BLen      Lengths of the branches, vector
% BAng      Branching angles of the branches, vector
% FCB       First cylinders of the branches, vector
%
% Additional outputs:
% TreeData  Vector containing basic tree attributes from the model
% BSeg      Segment of the branch, vector  (not every segment forms a branch)
% BChi      Child branches, cell-array
% CiB       Cylinders in the branches, cell-array
% CChi      Children cylinders of each cylinder, cell array
% CiS       Cylinders forming each segment, cell array
% Added     Logical vector indicating cylinders that are added to fill gaps
% P         Filtered point cloud, matrix
% Bal       Cover sets, cell array
% Cen       Center points of the cover sets, vector
% Nei       Neighboring cover sets, cell array
% Segs      Tree segments, cell array
% SPar      Parent segment of each segment, vector
% SChi      Child segments of each segment, cell array
% SDir      Direction lines of the segment bases, matrix
% SoC       Segments the cylinders belong, vector
% Base      Base of the tree
% Forb      Cover sets not part of the tree, vector

global RUNTIME;
RUNTIME = 1800;

% Names of the steps to display
name = ['Filtering  ';
        'Cover sets ';
        'Tree sets  ';
        'Segments   ';
        'Cylinders  '];
    
disp('---------------')
disp(string)
str = ['dmin0 = ',num2str(dmin0),', rcov0 = ',num2str(rcov0),', nmin0 = ',num2str(nmin0)];
disp(str)
str = ['dmin = ',num2str(dmin),', rcov = ',num2str(rcov),...
    ', nmin = ',num2str(nmin),', lcyl = ',num2str(lcyl)];
disp(str)
disp('Progress:')
tot = 0;

%%Read from .txt file and convert to .mat for use (added by Zane)
str = [P, '.txt'];
str2 = ['pointcloud/', P, '.mat'];
P = importdata(str);
save(str2, 'P');

%% Make the point cloud into proper form
% only 3-dimensional data
if size(P,2) > 3
    P = P(:,1:3);
end

% Only double precision data
C = class(P);
if C == 'single'
    P = double(P);
end

% If units in the point cloud are in millimeters,
% change them to meters
if max(P(:,3))-min(P(:,3)) > 1000
    P = P/1000;
    dmin = dmin/1000;
    rcov = rcov/1000;
end


%% Filtering
if nargin > 10
    tic
    str = ['FILTERING...', string];
    disp(str)
    %disp('FILTERING...')
    Pass = filtering(P,rfil1,nfil1,rfil2,nfil2);
    P = P(Pass,:);
    t = toc;
    tot = tot+t;
    [tmin,tsec] = sec2min(t);
    [Tmin,Tsec] = sec2min(tot);
    str = [name(1,:),' ',num2str(tmin),' min ',num2str(tsec),...
        ' sec',', total: ',num2str(Tmin),' min ',num2str(Tsec),' sec'];
    disp(str)
end


%% Cover sets
tic
disp('GENERATING COVER...')
[Bal,Nei,Cen] = cover_sets(P,dmin0,rcov0,nmin0);
t = toc;
tot = tot+t;
[tmin,tsec] = sec2min(t);
[Tmin,Tsec] = sec2min(tot);
str = [name(2,:),' ',num2str(tmin),' min ',num2str(tsec),...
    ' sec',', total: ',num2str(Tmin),' min ',num2str(Tsec),' sec'];
disp(str)


%% Tree sets
tic
disp('DETERMINING TREE SETS...')
[Base,Forb,Nei] = tree_sets(P,Cen,Bal,Nei,dmin0,NoGround);
t = toc;
tot = tot+t;
[tmin,tsec] = sec2min(t);
[Tmin,Tsec] = sec2min(tot);
str = [name(3,:),' ',num2str(tmin),' min ',num2str(tsec),...
    ' sec',', total: ',num2str(Tmin),' min ',num2str(Tsec),' sec'];
disp(str)


%% Segmenting
tic
disp('SEGMENTING...')
[Segs,SPar,SChi] = segments(P,Bal,Nei,Cen,Base,Forb,dmin0);

[Segs,SChi] = correct_segments(P,Cen,Segs,SPar,SChi,Bal,dmin0,1,1);
t = toc;
tot = tot+t;
[tmin,tsec] = sec2min(t);
[Tmin,Tsec] = sec2min(tot);
str = [name(4,:),' ',num2str(tmin),' min ',num2str(tsec),...
    ' sec',', total: ',num2str(Tmin),' min ',num2str(Tsec),' sec'];
disp(str)


%% New cover sets
tic
disp('GENERATING NEW COVER...')
% Determine relative size of new cover sets and use only tree points
[P,RS] = relative_size(P,Cen,Bal,Segs,SChi);
clearvars Bal Cen Segs SChi

[Bal,Nei,Cen] = cover_sets(P,dmin,rcov,nmin,RS);
t = toc;
tot = tot+t;
[tmin,tsec] = sec2min(t);
[Tmin,Tsec] = sec2min(tot);
str = [name(2,:),' ',num2str(tmin),' min ',num2str(tsec),...
    ' sec',', total: ',num2str(Tmin),' min ',num2str(Tsec),' sec'];
disp(str)


%% Tree sets
tic
disp('DETERMINING TREE SETS...')
[Base,Forb,Nei] = tree_sets(P,Cen,Bal,Nei,dmin,1);
t = toc;
tot = tot+t;
[tmin,tsec] = sec2min(t);
[Tmin,Tsec] = sec2min(tot);
str = [name(3,:),' ',num2str(tmin),' min ',num2str(tsec),...
    ' sec',', total: ',num2str(Tmin),' min ',num2str(Tsec),' sec'];
disp(str)


%% Segmenting
tic
disp('SEGMENTING...')
[Segs,SPar,SChi] = segments(P,Bal,Nei,Cen,Base,Forb,dmin);

[Segs,SChi,SPar] = correct_segments(P,Cen,Segs,SPar,SChi,Bal,dmin,1);
t = toc;
tot = tot+t;
[tmin,tsec] = sec2min(t);
[Tmin,Tsec] = sec2min(tot);
str = [name(4,:),' ',num2str(tmin),' min ',num2str(tsec),...
    ' sec',', total: ',num2str(Tmin),' min ',num2str(Tsec),' sec'];
disp(str)


%% Cylinders
tic
disp('CONSTRUCTING CYLINDER MODEL...')
[Rad,Len,Axe,Sta,CPar,CExt,CChi,SoC,CiS,Segs,SPar,SChi,Added] = ...
    cylinders(P,Bal,Cen,Segs,SPar,SChi,lcyl);

t = toc;
tot = tot+t;
[tmin,tsec] = sec2min(t);
[Tmin,Tsec] = sec2min(tot);
str = [name(5,:),' ',num2str(tmin),' min ',num2str(tsec),...
    ' sec',', total: ',num2str(Tmin),' min ',num2str(Tsec),' sec'];
disp(str)


tic
%% Determine the branches
[BoC,BOrd,BPar,BChi,CiB,BVol,BLen,BAng,FCB,BSeg,BHei,BAzi,BDia] = ...
    branches(Segs,SChi,CiS,CChi,CPar,Rad,Len,Axe,Sta,Added);


%% Compute and display model attributes
T = Segs{1};
T = vertcat(T{:});
T = vertcat(Bal{T});
trunk = P(T,:);
[TreeData,Taper,BODist,BSDist,CODist,CSDist,strX,strS,BO] = ...
    tree_data(Rad,Len,Sta,Axe,BOrd,CiB,BVol,BLen,trunk,string);
t = toc;
tot = tot+t;
[Tmin,Tsec] = sec2min(tot);

%save(string,'TreeData','BODist','BSDist','strX','BO','Taper',...
%   'CODist','CSDist','strS','string','Rad','Len','Sta','Axe','dmin','dmin0',...
%   'rcov','rcov0','nmin','nmin0','lcyl','Tmin','Tsec')

str = ['results/',string];

save(str,'TreeData','Sta','Axe','Rad','Len','CPar','CExt','BoC','BOrd','BPar', ...
     'BVol','BLen','BAng','BSeg','FCB','BChi','CiB','CChi','CiS','Added','P','Bal', ...
     'Cen','Nei','Segs','SPar','SChi','SoC','Base','Forb');


close all

%publish('make_report','pdf');

%str = ['results/results_report_',string,'.pdf'];
%movefile('report/html/make_report.pdf',str)


%% Plot models
% plot_branch_structure(P,Bal,Segs,BSeg,BChi,1,3,0,1);
% plot_cylinder_model(Rad,Len,Axe,Sta,2,20,0.3)
% plot_branch_structure(P,Bal,Segs,BSeg,BChi,1,8,0,1);
% 
% plot_tree_structure(P,Bal,Segs,SChi,1,2,0,1)
% hold on
% plot_cylinder_model(Rad,Len,Axe,Sta,1,20,1)
% hold off
 
%plot_branching_structure(BoC,Rad,Len,Axe,Sta,2,20,1)


%% Save the cylinder, branch, tree and distributions data in matlab and text-files
Rad = single(Rad);
Len = single(Len);
Sta = single(Sta);
Axe = single(Axe);
BVol = single(BVol);
BLen = single(BLen);
BAng = single(BAng);

rad = round(10000*Rad)/10000;
len = round(10000*Len)/10000;
sta = round(10000*Sta)/10000;
axe = round(10000*Axe)/10000;
bvol = round(1000*BVol)/1000;
blen = round(1000*BLen)/1000;
bang = round(1000*BAng)/1000;
CylData = single([rad len sta axe CPar CExt BoC Added]);
BranchData = single([BOrd BPar bvol blen bang BHei BAzi BDia]);

str = ['results/ModelData_',string];

save(str,'CylData','BranchData','TreeData');


% plot pointcloud and cylinder model ADDED BY KIM
% 
%grid off 
% plot_branching_structure(BoC,Rad,Len,Axe,Sta,1,20,0.8)
% hold on
% scatter3(P(:,1),P(:,2),P(:,3),0.001,'red')

%
%scatter3(P(:,1),P(:,2),P(:,3),0.001,P(:,3))
%axis equal
% 
%  hold off
% 

%plot_branching_structure(BoC,Rad,Len,Axe,Sta,2,20,1)

% this is the old code, doesn't do partial cylinders. Use one below
% Bot = min(Sta(:,3));
% Top = max(Sta(:,3)+Len.*Axe(:,3));
% H = ceil(Top-Bot);  % number of 1m height layers
% Vol = zeros(H,1);  % volumes of the bins
% for i = 1:H
% 	I = Sta(:,3) >= Bot+i-1;
% 	J = Sta(:,3) < Bot+i;
% 	I = I&J;
% 	Vol(i,1) = sum(pi*Len(I).*Rad(I).^2)*1000;
% end
% str = ['results/VolH_old',string,'.txt'];
% fid = fopen(str, 'wt');
% fprintf(fid, [repmat('%g\t', 1, size(TreeData,2)-1) '%g\n'], Vol.');
% fclose(fid);



%function vol_bins = heigth_bins(Rad,Len,Axe,Sta,bin_height,string)

% assumption: length of the cylinder shorter than bin_height
%bin_height=1
%Bot = min(Sta(:,3)); % bottom height
%Top = max(Sta(:,3)+Len.*Axe(:,3));  % top height
%H = ceil((Top-Bot)/bin_height);  % number of bin_height m layers
%vol_bins = zeros(H,1);  % volumes of the height bins
%n = length(Rad); % number of cylinders
%Sta(:,3) = Sta(:,3)-Bot; % move bottom height to zero
%for i = 1:H
%    vol = 0;
%    for j = 1:n
%        S = Sta(j,3);
%        T = S+Len(j)*Axe(j,3);
%        if T > S && S >= (i-1)*bin_height && S < i*bin_height % bottom inside the bin
%            t = (i*bin_height-S)/(Len(j)*Axe(j,3));
%            if t > 1 % top inside the bin
%                vol = vol+Rad(j)^2*Len(j);
%            else % top in the higher bin
%                vol = vol+t*Rad(j)^2*Len(j);
%            end
%        elseif T > S && T >= (i-1)*bin_height && T < i*bin_height % top inside the bin and bottom in the lower bin
%            t = (T-(i-1)*bin_height)/(Len(j)*Axe(j,3));
%            vol = vol+t*Rad(j)^2*Len(j);
%        elseif T <= S && T >= (i-1)*bin_height && T < i*bin_height % top inside the bin
%            t = (i*bin_height-T)/(Len(j)*abs(Axe(j,3)));
%            if t > 1 % bottom inside the bin
%                vol = vol+Rad(j)^2*Len(j);
%            else % bottom in the higher bin
%                vol = vol+t*Rad(j)^2*Len(j);
%            end
%        elseif T <= S && S >= (i-1)*bin_height && S < i*bin_height % bottom inside the bin and top in the lower bin
%            t = (S-(i-1)*bin_height)/(Len(j)*abs(Axe(j,3)));
%            vol = vol+t*Rad(j)^2*Len(j);
%        end
%    end
%    vol_bins(i) = 1000*pi*vol;
%end
%str = ['results/VolH',string,'.txt'];
%fid = fopen(str, 'wt');
%fprintf(fid, [repmat('%g\t', 1, size(vol_bins,2)-1) '%g\n'], vol_bins.');
%fclose(fid);





% 
% I = Rad > 0.05;
% LstemV=sum(pi*Len(I).*Rad(I).^2)*1000
% str = ['results/LstemV_data_',string,'.txt'];
% fid = fopen(str, 'wt');
% fprintf(fid,'%f%c',LstemV);
% fclose(fid);
% I = Rad > 0.025;
% UstemV=(sum(pi*Len(I).*Rad(I).^2)*1000)-LstemV
% str = ['results/UstemV_data_',string,'.txt'];
% fid = fopen(str, 'wt');
% fprintf(fid,'%f%c',UstemV);
% fclose(fid);
% I = Rad > 0;
% Can=(sum(pi*Len(I).*Rad(I).^2)*1000)-LstemV-UstemV
% str = ['results/Can_data_',string,'.txt'];
% fid = fopen(str, 'wt');
% fprintf(fid,'%f%c',Can);
% fclose(fid);

% matlab-format (.mat)
%str = ['results/cyl_data_',string];
%save(str,'CylData')

%str = ['results/branch_data_',string];
%save(str,'BranchData')

%str = ['results/tree_data_',string];
%save(str,'TreeData')

%str = ['results/distributions_',string];
%save(str,'Taper','BODist','BSDist','CODist','CSDist')

% text-format (.txt)
str = ['results/cyl_data_',string,'.txt'];
fid = fopen(str, 'wt');
fprintf(fid, [repmat('%g\t', 1, size(CylData,2)-1) '%g\n'], CylData.');
fclose(fid);

%str = ['results/branch_data_',string,'.txt'];
%fid = fopen(str, 'wt');
%fprintf(fid, [repmat('%g\t', 1, size(BranchData,2)-1) '%g\n'], BranchData.');
%fclose(fid);

%str = ['results/tree_data_',string,'.txt'];
%fid = fopen(str, 'wt');
%fprintf(fid, [repmat('%g\t', 1, size(TreeData,2)-1) '%g\n'], TreeData.');
%fclose(fid);
