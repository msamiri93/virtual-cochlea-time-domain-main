function AMat = make_AMatrix_FSI_CFLD(Nd,FLD,AMat,MP)

    CFLD = FLD.CFLD;
    ndof = 6; 

    % force on the structue from Corti pressure

    pdof = ((1:CFLD.Nd.N)-1)*CFLD.Nd.NV + 3; %pressure dof of Corti
 
    [~,idx] = sort(CFLD.Nd.x(CFLD.Nd.ind_Radi,1),'ascend');
    pdof_FSI = pdof(CFLD.Nd.ind_Radi(idx));
    partsAuc = zeros((size(CFLD.nBot,2)+size(CFLD.nTop,2))*1*4*(length(FLD.nBot)-1),4); % [2x2] matrices for 2 surfaces with 1 interacting node
    cntAuc = 0;
%     pAreaB = MP.dZ*2*Nd.X(CFLD.nBot(indbc_p,1))*CFLD.dcf_23; % factor of 2 is for getting the width of BM from the half point coordiante; area of effect for the Corti fluid should be reconsidered
%     pAreaR = MP.dZ*2*Nd.X(CFLD.nTop(indbc_p,1))*CFLD.dcf_23;

    pArea = MP.dZ*sqrt(4*CFLD.Nd.A(idx)/pi);
    % Refer to Pub 1 Supporting Materials.docx Eq.(10) (11) and (12)

    fdir = 2; % dof of load on structure

    nn = size(CFLD.nBot,2);
    for in = 1:nn
        % BM
        krow = (CFLD.nBot(:,in)-1)*ndof + fdir;
        kcol = pdof_FSI; % pressure dof from each node

        [partsAuc,cntAuc] = fillT(krow(2:end),kcol(1:end-1),-pArea(2:end)*MP.key_FSI/6/nn,partsAuc,cntAuc);
        [partsAuc,cntAuc] = fillT(krow(2:end),kcol(2:end),-pArea(2:end)*MP.key_FSI/3/nn,partsAuc,cntAuc);
        [partsAuc,cntAuc] = fillT(krow(1:end-1),kcol(2:end),-pArea(1:end-1)*MP.key_FSI/6/nn,partsAuc,cntAuc);
        [partsAuc,cntAuc] = fillT(krow(1:end-1),kcol(1:end-1),-pArea(1:end-1)*MP.key_FSI/3/nn,partsAuc,cntAuc);

    end
    
    nn = size(CFLD.nTop,2);
    for in = 1:nn
        krow = (CFLD.nTop(:,in)-1)*ndof + fdir; % column index
        kcol = pdof_FSI; % pressure dof from each node

        [partsAuc,cntAuc] = fillT(krow(2:end),kcol(1:end-1),pArea(2:end)*MP.key_FSI/6/nn,partsAuc,cntAuc);
        [partsAuc,cntAuc] = fillT(krow(2:end),kcol(2:end),pArea(2:end)*MP.key_FSI/3/nn,partsAuc,cntAuc);
        [partsAuc,cntAuc] = fillT(krow(1:end-1),kcol(2:end),pArea(1:end-1)*MP.key_FSI/6/nn,partsAuc,cntAuc);
        [partsAuc,cntAuc] = fillT(krow(1:end-1),kcol(1:end-1),pArea(1:end-1)*MP.key_FSI/3/nn,partsAuc,cntAuc);
   
    end

    Auc = sparse(partsAuc(:,1),partsAuc(:,2),partsAuc(:,3),Nd.tdof,CFLD.tdof);
    rdof = Nd.BC;
    cdof = CFLD.BC;
    AMat.Auc = Auc(rdof,cdof);

%     Qcc = CFLD.Q(CFLD.BC,CFLD.BC);
%     AMat.Qcc_O = Qcc;  
end