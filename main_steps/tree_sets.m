function [Base,Forb,Nei] = tree_sets(P,Cen,Bal,Nei,dmin,NoGround)

% ---------------------------------------------------------------------
% TREE_SETS.M       Determines the base of the trunk and the cover sets 
%                   belonging to the tree, updates the neighbor-relation
%
% Version 1.3
% Author        Pasi Raumonen
% Created       21 February 2014
%
% This software cannot be distributed without the permission of P. Raumonen
% and it is only for non-commercial use. These restrictions holds also for 
% the modules or subprograms of the software:
% CONNECTED_COMPONENTS.M, PARTITION_OF_POINT_CLOUD.M.
% ---------------------------------------------------------------------

% Determines the cover sets that belong to the tree. Determines also the
% base of the tree and updates the neighbor-relation such that all of the
% tree is connected, i.e., the cover sets belonging to the tree form a
% single connected component.

% Inputs:
% P         Point cloud
% Cen       Center points of the cover sets, (n_sets x 1)-vector
% Bale      Extended cover sets, (n_sets x 1)-cell
% Nei       Neighboring cover sets, (n_sets x 1)-cell
% Dir       Direction vectors of the cover sets (n_sets x 3)-maxtrix
% Dim       Dimensionality values of the sets, (n_sets x 3)-matrix
% dmin      Minimum diameter of the cover sets
%
% Outputs:
% Base      Base of the trunk (vector)
% Forn      Cover sets not part of the tree
% Bale      Updated extended cover sets
% Nei       Updated neigbors


nb = max(size(Bal));   % number of cover sets
Fal = false(nb,1);
Ind = (1:1:nb)';
Ce = P(Cen,1:3); % Coordinates of the center points
Hmin = min(Ce(:,3));
Height = max(Ce(:,3))-Hmin;

if NoGround
    % No ground in the point cloud
    % The base is the lowest part
    H = min([0.1 0.02*Height]);
    I = Ce(:,3) < Hmin+H;
    Base = Ind(I);
    Forb = false(nb,1);
    
    % Determine the Trunk
    Trunk = Fal;
    Trunk(Base) = true;
    Added = Trunk;
    Added(vertcat(Nei{Added})) = true;
    Added(Trunk) = false;
    Trunk(Added) = true;
    L = 0.05*Height;
    H = max(Ce(Trunk,3))-L;
    while any(Added)
        Added(vertcat(Nei{Added})) = true;
        Added(Trunk) = false;
        A = Ind(Added);
        I = Ce(Added,3) < H;
        Added(A(I)) = false;
        Trunk(Added) = true;
        H = max(Ce(A(~I),3))-L;
    end
else
    % Determine the base from the "height" and "density" of cover sets
    % by projecting the sets to the xy-plane
    
    % The vertices of the rectangle containing C
    Min = double(min(Ce));
    Max = double(max(Ce(:,1:2)));
    
    % Number of rectangles with edge length "E" in the plane
    E = 0.015*Height;
    n = double(ceil((Max(1:2)-Min(1:2))/E)+1);
    
    % Calculates the rectangular-coordinates of the points
    px = floor((Ce(:,1)-Min(1))/E)+1;
    py = floor((Ce(:,2)-Min(2))/E)+1;
    
    % Sorts the points according a lexicographical order
    S = [px py-1]*[1 n(1)]';
    [S,I] = sort(S);
    
    Partition = cell(n(1),n(2));
    h = zeros(n(1),n(2)); % "height" of the cover sets in the squares
    d = h;  % density of the cover sets in the squares
    b = h;
    p = 1; % The index of the point under comparison
    while p <= nb
        t = 1;
        while (p+t <= nb) && (S(p) == S(p+t))
            t = t+1;
        end
        q = I(p);
        J = I(p:p+t-1);
        Partition{px(q),py(q)} = J;
        p = p+t;
        K = ceil(10*(Ce(J,3)-Min(3)+0.01)/(Height-0.01));
        B = K <= 2;
        K = unique(K);
        h(px(q),py(q)) = length(K)/10;
        d(px(q),py(q)) = t;
        b(px(q),py(q)) = nnz(B);
    end
    d = d/max(max(d));  % normalize
    b = b/max(max(b));
    
    f = d.*h.*b;  % function whose maximum determines location of the trunk
    % smooth the function by averaging over 8-neighbors
    for i = 2:n(1)-1
        for j = 2:n(2)-1
            f(i,j) = mean(mean(f(i-1:i+1,j-1:j+1)));
        end
    end
    f = f/max(max(f));
    
    % Trunk location is around the maximum f-value 
    I = f > 0.5;
    Tru = Partition(I); % squares that contain the trunk
    Tru = vertcat(Tru{:});
    h = min(Ce(Tru,3));
    I = Ce(Tru,3) > h+0.025*Height;
    J = Ce(Tru,3) < h+0.15*Height;
    I = I&J; 
    Tru = Tru(I);  % "slice" between 0.5-2m that should contain trunk 
    Tru = union(Tru,vertcat(Nei{Tru})); % Expand
    Tru = union(Tru,vertcat(Nei{Tru})); % Expand
    
    % Define connected components of Tru and select the largest component
    [Comp,CS] = connected_components(Nei,Tru,10,Fal);
    [~,I] = max(CS);
    Tru = Comp{I};
    
    % Fit cylinder to Tru
    I = Ce(Tru,3) < h+0.15*Height;
    Tru = Tru(I);
    Tru = union(Tru,vertcat(Nei{Tru}));
    Points = Ce(Tru,:);
    AP0 = mean(Points)';
    CA0 = [0 0 1]';
    R0 = mean(distances_to_line(Points,CA0,AP0));
    [AP, CA, R, d] = lscylinder(Points, AP0, CA0, R0, 0.1, 0.1);
    
    % Remove far away points and fit new cylinder
    I = d < 0.1*R;
    Points = Points(I,:);
    Tru = Tru(I);
    [AP, CA, R] = lscylinder(Points, AP, CA, R, 0.1, 0.1);
    H = Points*CA;
    hmin = min(H);
    hpoint = CA'*AP;
    AP = AP-(hpoint-hmin)*CA;
    
    % Expand downwards as long as the radius does not grow too much
    q = Tru;
    r = R;
    again = true;
    while again
        I = Ce(Tru,3) < h+0.05*Height;
        Tru = Tru(I);
        Tru = union(Tru,vertcat(Nei{Tru}));
        Points = Ce(Tru,:);
        AP0 = AP;
        CA0 = CA;
        R0 = R;
        [AP, CA, R, d] = lscylinder(Points, AP0, CA0, R0, 0.1, 0.1);
        I = d < 1.1*R;
        Tru = Tru(I);
        if (abs(R0-R) < 1e-5) || (R > 1.2*r && min(Points(:,3)) < h+0.01*Height)
            again = false;
        end
    end
    
    % Compute normals to the Tru sets and then remove those that are too
    % parallel to the cylinder axis
    nt = length(Tru);
    I = false(nt,1);
    for i = 1:nt
        B = vertcat(Bal{[i; Nei{Tru(i)}]});
        if length(B) > 3
            X = cov(P(B,:));
            [U,~,~] = svd(X);
            if abs(U(:,3)'*CA) < 0.5
                I(i) = true;
            end
        else
            I(i) = true;
        end
    end
    Tru = Tru(I);
    
    % Determine the base
    Bot = min(Ce(Tru,3));
    J = Ce(Tru,3) < Bot+0.075*Height;%0.025*Height;
    Base = Tru(J);
    Tru = union(Tru,q);
    bottom = min(Ce(Base,3));
    
    % Determine "trunk" by going up
    Trunk = Fal;
    Trunk(Tru) = true;
    G = setdiff(vertcat(Nei{Base}),Tru);
    G = setdiff(union(G,vertcat(Nei{G})),Tru);
    Added = Trunk;
    Added(vertcat(Nei{Added})) = true;
    Added(Trunk) = false;
    Added(G) = false;
    Trunk(Added) = true;
    L = 0.05*Height;
    H = max(Ce(Trunk,3))-L;
    while any(Added)
        Added(vertcat(Nei{Added})) = true;
        Added(Trunk) = false;
        A = Ind(Added);
        I = Ce(Added,3) < H;
        Added(A(I)) = false;
        Trunk(Added) = true;
        H = max(Ce(A(~I),3))-L;
    end
    
    % Determine ground, i.e. the "Forb", by expanding as much as possible
    Forb = Fal;
    Forb(G) = true;
    Forb(Base) = false;
    Added = Forb;
    while any(Added)
        Added(vertcat(Nei{Added})) = true;
        Added(Forb) = false;
        Added(Trunk) = false;
        Forb(Added) = true;
    end
    
    % Try to expand the "Forb" more by adding all the bottom sets
    G = Ce(:,3) < Bot+0.05*Height;
    Forb(G) = true;
    Forb(Trunk) = false;
end

%% Update the neighbor-relation near the Trunk
% Do the update near the trunk twise
for k = 1:2
    if k == 2
        %% Expand Trunk as much as possible
        Trunk(Forb) = false;
        Trunk(Base) = false;
        Added = Trunk;
        while any(Added)
            Added(vertcat(Nei{Added})) = true;
            Added(Trunk) = false;
            Added(Forb) = false;
            Added(Base) = false;
            Trunk(Added) = true;
        end
    end
    
    % Select only the cover sets nearest to the "Trunk"
    T = Ind(Trunk);
    n = length(T);
    EdgeLength = min(0.15,5*dmin);
    [Par,CC,info] = partition_of_point_cloud(Ce,EdgeLength);
    
    I = false(info(4),info(5),info(6)); % cover sets near the "Trunk"
    for j = 1:n
        J = T(j);
        I(CC(J,1)-1:CC(J,1)+1,CC(J,2)-1:CC(J,2)+1,CC(J,3)-1:CC(J,3)+1) = true;
    end
    N = Fal;
    N(vertcat(Par{I})) = true;
    N(T) = false;   % Only the non-trunk cover sets
    N = Ind(N);
    
    if ~isempty(N)
        % Determine separate components of "N", cover sets near the trunk
        Comps = connected_components(Nei,N,1,Fal);
        nc = size(Comps,1);
        if ~NoGround
            % use only sets high enough
            I = true(nc,1);
            for i = 1:nc
                C = Comps{i};
                if min(Ce(C,3)) < bottom+0.3*Height
                    I(i) = false;
                end
            end
            Comps = Comps(I);
            nc = nnz(I);
        end
        
        % Expand each component
        for i = 1:nc
            C = Comps{i};
            C = union(C,vertcat(Nei{C}));
            C = union(C,vertcat(Nei{C}));
            C = union(C,vertcat(Nei{C}));
            C = union(C,vertcat(Nei{C}));
            C = union(C,vertcat(Nei{C}));
            Comps{i} = C;
        end
        
        N = Fal;
        N(vertcat(Comps{:})) = true;
        N(T) = false;    % Only the non-trunk cover sets
        N = Ind(N);
        
        if ~isempty(N)
            % Determine separate components of "N", cover sets near the trunk
            [Comps,CompSize] = connected_components(Nei,N,1,Fal);
            nc = length(CompSize);
            I = true(nc,1);
            if ~NoGround
                % use only sets high enough
                for i = 1:nc
                    C = Comps{i};
                    if min(Ce(C,3)) < bottom+0.3*Height
                        I(i) = false;
                    end
                end
                Comps = Comps(I);
                nc = nnz(I);
            end
            
            % Check the components and possibly update the neighbor-relation
            K = false(nc,1);
            for i = 1:nc
                C = Comps{i};
                J = Trunk(C);
                L = Trunk(vertcat(Nei{C}));
                
                if any(J) || any(L)
                    %disp('yhteydessa')
                else
                    NC = length(C);
                    % Select only the cover sets the nearest to the component
                    I = false(info(4),info(5),info(6));
                    for j = 1:NC
                        J = C(j);
                        I(CC(J,1)-1:CC(J,1)+1,CC(J,2)-1:CC(J,2)+1,CC(J,3)-1:CC(J,3)+1) = true;
                    end
                    
                    M = Fal;
                    M(vertcat(Par{I})) = true;
                    M(~Trunk) = false; % The nearest "Trunk" cover sets
                    M = Ind(M);
                    
                    d = pdist2(Ce(C,:),Ce(M,:));
                    if NC == 1 && length(M) == 1
                        dt = d;
                        I = C;
                        J = M;
                    elseif NC == 1
                        [dt,J] = min(d);
                        I = C;
                        J = M(J);
                    elseif length(M) == 1
                        [dt,I] = min(d);
                        J = M;
                        I = C(I);
                    else
                        [d,I] = min(d);
                        [dt,J] = min(d);
                        I = I(J);
                        I = C(I);
                        J = M(J);
                    end
                    
                    if (NC > 5) || ((NC <= 5) && (dt < 3*dmin))
                        
                        K(i) = true;
                        Nei{I} = [Nei{I}; J];
                        Nei{J} = [Nei{J}; I];
                    end
                    
                end
            end
        end
    end
end

%% Update the neighbor relation for whole tree
% Define the expanded trunk or "Tree" again
Tree = Fal;
Tree(vertcat(Nei{Base})) = true;
Tree(Forb) = false;
Tree(Base) = false;
Added = Tree;
while any(Added)
    Added(vertcat(Nei{Added})) = true;
    Added(Tree) = false;
    Added(Forb) = false;
    Added(Base) = false;
    Tree(Added) = true;
end

Other = ~Fal;  % sets not yet connected to "Tree"
Tree(Base) = false;
Other(Tree) = false;
Other(Forb) = false;
k0 = min(10,ceil(0.4/dmin));  % multiplier for dmin to determine the closeness of components
k = k0;
Cmin = ceil(0.2/dmin);  % minimum accepted component size
Cmin = 0;
W = 1000;  % If a component has more than W sets, then divide the distance calculations into smaller parts
while any(Other)
    npre = nnz(Other);
    again = true;
    
    % Partition the centers of the cover sets into cubes
    [Par,CC,info] = partition_of_point_cloud(Ce,k*dmin);
    
    while any(Other) && again
        
        % Determine the components of "Other", not connected to "Tree"
        Comps = connected_components(Nei,Other,1,Fal);
        nc = size(Comps,1);
        
        % Label the cover sets by their components
        CL = zeros(nb,1,'int32');
        CL(Tree) = nc+1;
        CL(Forb) = nc+3;
        for i = 1:nc
            CL(Comps{i}) = i;
        end
        
        % Try to connect components to "Tree", "Other", or "Forb"
        for i = 1:nc
            C = Comps{i};
            NC = length(C);
            
            % Select only the cover sets the nearest to the component
            I = false(info(4),info(5),info(6));
            for j = 1:NC
                I(CC(C(j),1)-1:CC(C(j),1)+1,CC(C(j),2)-1:CC(C(j),2)+1,CC(C(j),3)-1:CC(C(j),3)+1) = true;
            end
            N = Fal;
            N(vertcat(Par{I})) = true;
            N(C) = false;  % The nearest cover sets
            N = Ind(N);

            L = CL(N);  % The component labels of the nearest cover sets
            Tre = L == nc+1;  % "Tree" sets
            F = L == nc+3;  % "Forb" sets
            
            O = N(~(Tre|F));
            Tre = N(Tre);
            F = N(F);
            
            % Determine the closest sets for "Other", "Tree" and "Forb"
            if ~isempty(O)
                if NC > W
                    m = ceil(NC/W);
                    Dmin = 100;
                    for j = 1:m
                        if j < m
                            d = pdist2(Ce(C((j-1)*W+1:j*W),:),Ce(O,:));
                            D = min(min(d));
                            if D < Dmin
                                Dmin = D;
                                dis = d;
                            end
                        else
                            d = pdist2(Ce(C((j-1)*W+1:end),:),Ce(O,:));
                            D = min(min(d));
                            if D < Dmin
                                Dmin = D;
                                dis = d;
                            end
                        end
                    end
                    d = dis;
                else
                    d = pdist2(Ce(C,:),Ce(O,:));
                end
                if NC == 1 && length(O) == 1
                    do = d;
                    IO = 1;
                    JO = 1;
                elseif NC == 1
                    [do,JO] = min(d);
                    IO = 1;
                elseif length(O) == 1
                    [do,IO] = min(d);
                    JO = 1;
                else
                    [d,IO] = min(d);
                    [do,JO] = min(d);
                    IO = IO(JO);
                end
            else
                do = 9;
            end
            
            if ~isempty(Tre)
                if NC > W
                    m = ceil(NC/W);
                    Dmin = 100;
                    for j = 1:m
                        if j < m
                            d = pdist2(Ce(C((j-1)*W+1:j*W),:),Ce(Tre,:));
                            D = min(min(d));
                            if D < Dmin
                                Dmin = D;
                                dis = d;
                                CTe = C((j-1)*W+1:j*W);
                            end
                        else
                            d = pdist2(Ce(C((j-1)*W+1:end),:),Ce(Tre,:));
                            D = min(min(d));
                            if D < Dmin
                                Dmin = D;
                                dis = d;
                                CTe = C((j-1)*W+1:end);
                            end
                        end
                    end
                    d = dis;
                else
                    d = pdist2(Ce(C,:),Ce(Tre,:));
                end
                if NC == 1 && length(Tre) == 1
                    dt = d;
                    IT = 1;
                    JT = 1;
                elseif NC == 1
                    [dt,JT] = min(d);
                    IT = 1;
                elseif length(Tre) == 1
                    [dt,IT] = min(d);
                    JT = 1;
                else
                    [d,IT] = min(d);
                    [dt,JT] = min(d);
                    IT = IT(JT);
                end
            else
                dt = 8;
            end
            
            if ~isempty(F)
                if NC > W
                    m = ceil(NC/W);
                    Dmin = 100;
                    for j = 1:m
                        if j < m
                            d = pdist2(Ce(C((j-1)*W+1:j*W),:),Ce(F,:));
                            D = min(min(d));
                            if D < Dmin
                                Dmin = D;
                                dis = d;
                            end
                        else
                            d = pdist2(Ce(C((j-1)*W+1:end),:),Ce(F,:));
                            D = min(min(d));
                            if D < Dmin
                                Dmin = D;
                                dis = d;
                            end
                        end
                    end
                    d = dis;
                else
                    d = pdist2(Ce(C,:),Ce(F,:));
                end
                df = min(d);
                if length(df) > 1
                    df = min(df);
                end
            else
                df = 10;
            end
            
            d = min([do dt]);
            
            % Determine what to do with the component
            if NC < Cmin && d > 0.3
                % Remove small isolated component
                Forb(C) = true;
                Other(C) = false;
                CL(C) = nc+3;
            elseif df <= do && df <= dt
                % Join the component to "Forb"
                Forb(C) = true;
                Other(C) = false;
                CL(C) = nc+3;
            elseif df == 10 && dt == 8 && do == 9
                % Isolated component, do nothing
            else
                if (dt <= do) || (dt <= 2*dmin)
                    % Join to "Tree"
                    if NC > W
                        C = CTe;
                    end
                    I = C(IT);
                    J = Tre(JT);
                    C = Comps{i};
                    Other(C) = false;
                    Tree(C) = true;
                    CL(C) = nc+2;
                else
                    % Join to other components
                    I = C(IO);
                    J = O(JO);
                end
                Nei{I} = [Nei{I}; J];
                Nei{J} = [Nei{J}; I];
            end
        end
        
        % If "Other" has decreased, do another check with same distance
        if nnz(Other) < npre
            again = true;
            npre = nnz(Other);
        else
            again = false;
        end
    end
    k = k+k0;
    Cmin = 3*Cmin;
end
if nnz(Other) > 0
    str = ['Problem!! ',numstr(nnz(Other)),' sets not included'];
    disp(str)
end
Forb(Base) = false;

%plot_segments(P,Bal,1,Base,Ind(Forb),Ind)