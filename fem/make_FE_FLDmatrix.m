function [AMat, FLD] = make_FE_FLDmatrix(MP,AMat,FLD,opt_cfld)
    % Input parameters
    Nd = FLD.Nd;
    El = FLD.El;

%     nnel = size(t,2); % number of nodes per element
    ned = 2; % number of nodes per element edge
    NN = El.type;
    ndof = 1; % number of dofs per node
    sdof = Nd.N*ndof; % total dofs
    
    partsMA = zeros(NN*NN*El.N,3);
    partsMQ = zeros(ned*ned*El.N,3);
    cntA = 0;
    cntQ = 0;
 
    bcBM = find(strncmp(FLD.name,'BM',2));
    [~,idx] = sort(FLD.Nd.x(bcBM,1),'ascend');
    bcBM = bcBM(idx);

    if eq(El.type,4)
        [~, ~, w] = set_Gauss_local_variables(2,2,El.NQ);
    else
        w = [];
    end
    for ie = 1:El.N

        knodes = FLD.El.node(ie,:);   % nodes of element kk
        A_e = D_matrix(ie,Nd,El,w);
        [partsMA, cntA] = fillM(ndof,knodes', A_e, partsMA, cntA);
    
        if opt_cfld == 1
            wrap = [knodes,knodes(1)];
            for id = 1:NN
                j1 = wrap(id); j2 = wrap(id+1);
                i1 = find(bcBM == j1,1);
                i2 = find(bcBM == j2,1);
                if (~isempty(i1) && ~isempty(i2))
                    xe1 = Nd.x(j1,1); ye1 = Nd.x(j1,2);
                    xe2 = Nd.x(j2,1); ye2 = Nd.x(j2,2);
                    edge = sqrt((xe2-xe1)^2+(ye2-ye1)^2);
                    kk = perm_scale(MP,mean(FLD.CFLD.Nd.Z([i1,i2])));
                    width = mean(FLD.CFLD.Nd.width([i1,i2]));
                    eqm = FLD.rho*kk*width*eqm2(edge); % (-1 normal unit inward)*( - pressure node on bottom of BM)
                    [partsMQ, cntQ] = fillM(ndof, [j1;j2], eqm, partsMQ, cntQ);
                end
            end
        end
    end

    if opt_cfld == 1
        partsMQ = partsMQ(any(partsMQ(:,[1,2]),2),:);    
        Q = sparse(partsMQ(:,1), partsMQ(:,2), partsMQ(:,3), sdof, sdof);
        AMat.Qpp_O = Q(FLD.BC,FLD.BC);
    end

    L = sparse(partsMA(:,1), partsMA(:,2), partsMA(:,3), sdof, sdof);
    FLD.L = L;
%     a = attenuation_function(MP,Nd.x(:,1),1,FLD.tdof,FLD.BC,1,3,10);
    AMat.App = L(FLD.BC,FLD.BC);
    
end

function eqm = eqm2(ds)

    % page 303 Pozrikidis, intro to finite and spectral element methods using matlab
    eqm(1,1) = ds/3;
    eqm(1,2) = ds/6;
    eqm(2,1) = ds/6;
    eqm(2,2) = ds/3;
    
end

function A_e = D_matrix(ie,Nd,El,w)
    
    A_e = zeros(El.type);    

    if eq(El.type,4)

%         global GaussQuad;        
%         w = GaussQuad.w;
        
        hs = El.elm_hs(:,:,ie);
    
        for iq = 1:El.NQ
            cf = hs(iq)*w(iq);
            gpsix = El.elm_gpsi(1:El.type,iq,ie);
            gpsiy = El.elm_gpsi(1 + El.type:end,iq,ie);
            A_e = A_e + cf*(gpsix*(gpsix.') + gpsiy*(gpsiy.'));
        end
    elseif eq(El.type,3)
        knodes = El.node(ie,:);   % nodes of element kk
        x = Nd.x(knodes,1); y = Nd.x(knodes,2);
        A = 1/2*det([1,x(1),y(1);1,x(2),y(2);1,x(3),y(3)]);
        A_e(1,1) = 1/4/A*((x(3)-x(2))^2+(y(2)-y(3))^2);
        A_e(1,2) = 1/4/A*((x(3)-x(2))*(x(1)-x(3)) + (y(2)-y(3))*(y(3)-y(1)));
        A_e(1,3) = 1/4/A*((x(3)-x(2))*(x(2)-x(1)) + (y(2)-y(3))*(y(1)-y(2)));
        A_e(2,2) = 1/4/A*((x(1)-x(3))^2 + (y(3)-y(1))^2);
        A_e(2,3) = 1/4/A*((x(1)-x(3))*(x(2)-x(1)) + (y(3)-y(1))*(y(1)-y(2)));
        A_e(3,3) = 1/4/A*((x(2)-x(1))^2 + (y(1)-y(2))^2);
        A_e(2,1) = A_e(1,2);
        A_e(3,1) = A_e(1,3);
        A_e(3,2) = A_e(2,3);
    end

end
