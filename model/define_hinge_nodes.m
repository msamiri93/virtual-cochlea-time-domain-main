function [Nd,El] = define_hinge_nodes(Nd,El)

    El.hinge = false(El.N,2);

    eOHB = strcmp('OHB',El.name); % OHC bundle and linear spring Nd1:bb, Nd2:e1
    El.hinge(eOHB,2) = 1;

   % %#% Hinge at the IPC root
   % eH = strcmp(Nd.name(El.Nd1),'CI') & strcmp(Nd.name(El.Nd2),'CI1');
   % El.hinge(eH,1) = 1;

   % #% Hinge at the OPC root
   eH = strcmp(Nd.name(El.Nd1),'AP') & strcmp(El.name, 'BMp');
   El.hinge(eH,1) = 1;

    %#% Hinge at the DC root
    eH = strcmp(Nd.name(El.Nd1),'AD') & strcmp(El.name, 'DCb'); 
    El.hinge(eH,1) = 1;

    dup_Nd = zeros(size(Nd.name));
    for ie = 1:El.N
        if sum(El.hinge(ie,:))>0
            nds = [El.Nd1(ie),El.Nd2(ie)];
            idx = find(El.hinge(ie,:));
            dup_Nd(nds(idx)) = nds(idx);
        end
    end
    idx = any(dup_Nd,1);
    h_name = cellfun(@(c)[c,'h'],Nd.name(:,idx),'UniformOutput',false);
    
    Nd.X = [Nd.X,Nd.X(:,idx)];
    Nd.Y = [Nd.Y,Nd.Y(:,idx)];
    Nd.Z = [Nd.Z,Nd.Z(:,idx)];
    Nd.name = [Nd.name,h_name];

    Nd = initialize_nodes(Nd);

    nd_con = zeros(sum(El.hinge(:)),2);

    icount = 1;
    for ie = 1:El.N
        if sum(El.hinge(ie,:))>0
            nds = [El.Nd1(ie),El.Nd2(ie)];
            idx = find(El.hinge(ie,:));
            r = mod(nds(idx),Nd.Nr); if r == 0, r = Nd.Nr; end
            c = ceil(nds(idx)/Nd.Nr);
            ch = find(strcmp(Nd.name(r,:),[Nd.name{r,c},'h']));
            nd_con(icount,:) = [nds(idx), r + Nd.Nr*(ch-1)];
            icount = icount + 1;
            if idx == 1
                El.Nd1(ie) = r + Nd.Nr*(ch-1);
            else
                El.Nd2(ie) = r + Nd.Nr*(ch-1);
            end
            if strncmp(El.name(ie),'DCb',3)
                El.Nd1(end+1) = nds(1);
                El.Nd2(end+1) = nds(2);
                El.name{end+1} = 'DCbc';
                El.hinge(end+1,:) = false(1,2);
            end
        end
    end
    El.N = numel(El.name);
    Nd.Master_Slave = nd_con;
end

function Nd = initialize_nodes(Nd)
        
    Nd.N = length(Nd.X);
    Nd.Nc = size(Nd.X,2);
    Nd.Nr = size(Nd.X,1);
        
    Nd.dx = zeros(size(Nd.X));
    Nd.dy = zeros(size(Nd.X));
    Nd.dz = zeros(size(Nd.X));
    Nd.rx = zeros(size(Nd.X));
    Nd.ry = zeros(size(Nd.X));
    Nd.rz = zeros(size(Nd.X));
    
    NDOF = 6;
    Nd.BC = ones(Nd.Nr* Nd.Nc, NDOF);
    Nd.N = Nd.Nr*Nd.Nc;
end