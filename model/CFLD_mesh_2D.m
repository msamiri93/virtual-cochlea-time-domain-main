function CFLD = CFLD_mesh_2D(MP,sNd,CFLD)

    xx = sNd.Z(CFLD.nBot);
    Nd.width = CFLD.Nd.width;
    Nd.Z = xx;
    Nd.nz = length(xx);


    % idx = [find(strcmp(sNd.name,'BB')),find(strcmp(sNd.name,'DD')),...
    %        find(strcmp(sNd.name,'AM')),find(strcmp(sNd.name,'A5')),...
    %        find(strcmp(sNd.name,'AP')),find(strcmp(sNd.name,'A3')),...
    %        find(strcmp(sNd.name,'A2')),find(strcmp(sNd.name,'A1')),...
    %        find(strcmp(sNd.name,'A0')),find(strcmp(sNd.name,'C1')),...
    %        find(strcmp(sNd.name,'C2')),find(strcmp(sNd.name,'C3')),...
    %        find(strcmp(sNd.name,'CC')),find(strcmp(sNd.name,'B1'))]; % entire OoC
    idx = [find(strcmp(sNd.name,'BB')),find(strcmp(sNd.name,'DD')),...
           find(strcmp(sNd.name,'AD')),find(strcmp(sNd.name,'AF')),...
           find(strcmp(sNd.name,'AP')),find(strcmp(sNd.name,'A1')),...
           find(strcmp(sNd.name,'A0')),find(strcmp(sNd.name,'CI')),...
           find(strcmp(sNd.name,'CI1')),find(strcmp(sNd.name,'CC')),...
           find(strcmp(sNd.name,'B1'))]; % entire OoC


    idx = idx(CFLD.idx,:);

    Nd.p_num = size(idx,2);
    Nd.perim = idx;

    A = zeros(Nd.nz,1);
    H = zeros(Nd.nz,2*Nd.p_num);

    for ni = 1:Nd.nz
        idx = Nd.perim(ni,:).';
        xx = [sNd.X(idx),sNd.Y(idx)];
            % formulation for approximate calculation of area change based on the shoelace formula
        A(ni,1) = polyarea(sNd.X(idx),sNd.Y(idx));
        H(ni,:) = -1*polyarea_change(xx); % the -1 is for the counterclock-wise direction of the node numbering, it could be reversed for convenience    
    end

    if MP.vMC == 1
        load('./assets/Corti_MC.mat','msh');
    else
        load('./assets/Corti_2D.mat','msh');
    end
    Nd.x = msh.POS(:,[1,2]);    
    El.node = msh.QUADS8(:,1:8);
    El.type = 8;
 
    El.N = size(El.node,1);
    Nd.N = size(Nd.x,1);
    Nd.H = H;
    Nd.A = A;
    Nd.radi = sqrt(A/pi);
    Nd.var_name = {'u','v','p','w'};
    Nd.NV = 4;

    eps = 1e-3;
    Nd.ind_Base = Nd.x(:,1) < min(Nd.x(:,1)) + eps;
    Nd.ind_Apex = Nd.x(:,1) > max(Nd.x(:,1)) - eps;
    Nd.ind_Symm = Nd.x(:,2) < min(Nd.x(:,2)) + eps;

    
    if MP.vMC == 1
        x = linspace(0,max(Nd.x(:,1)),401); % the number matches the # of nodes on the boundary
        y = interp1([0,max(Nd.x(:,1))],[25,25],x);
    else
        x = linspace(0,max(Nd.x(:,1)),1201);
        y = interp1([0,max(Nd.x(:,1))],[28.134,62.721],x);        
    end
    Nd.ind_Radi = knnsearch(Nd.x,[x(:),y(:)]);
    
    CFLD.Nd = Nd;
    CFLD.El = El;

end