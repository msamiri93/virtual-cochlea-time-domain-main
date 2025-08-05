function Aeo = make_Aeo(MP, OHC, Vnd)

    nv = 5;
    S = [OHC.S];
    OHC_Gs_max = [OHC.count].*[S.Gmax];

    odof = MP.dof(5).n;

    edof = MP.dof(3).n;

    % circuit model has 5 nodes per section
    rdof_node1 = 1:5:edof; % V1
    rdof_node2 = 2:5:edof; % V2

    Vnd = reshape(Vnd,nv,odof);
    V1 = Vnd(1,:);
    V2 = Vnd(2,:);

    col = 1:odof;

    Gsmax_by_V0 = OHC_Gs_max.*(V1 - V2);

    partsAeo = zeros(2*odof,4);
    cntAeo = 0;

    % equations 1 and 2 of the circuit containing stereociliary conductance
    [partsAeo,cntAeo] = fillT(rdof_node1,col,Gsmax_by_V0,partsAeo,cntAeo);
    [partsAeo,~] = fillT(rdof_node2,col,-Gsmax_by_V0,partsAeo,cntAeo);   

    Aeo = sparse(partsAeo(:,1),partsAeo(:,2),partsAeo(:,3),edof,odof);
end