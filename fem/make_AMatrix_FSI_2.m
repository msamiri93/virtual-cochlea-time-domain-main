function [AMat,FLD] = make_AMatrix_FSI_2(AMat,Nd,El,FLD,MP)

opt_cfld = MP.opt_cfld;
FLD.key_FSI_2Chambers = 3;

if opt_cfld == 1
    opt_mean_pressure = 0;
else
    opt_mean_pressure = 1;
end

ndof = 6;
bcTM = FLD.indTM; % MSA: TM is old naming convention, the RL is being used as an FSI boundary now!
[~,idx] = sort(FLD.Nd.x(bcTM,1),'ascend');
bcTM = bcTM(idx);
bcBM = FLD.indBM;
[~,idx] = sort(FLD.Nd.x(bcBM,1),'ascend');
bcBM = bcBM(idx);
surface = FLD.nBot;

fdir = 2; 

partsApu_O = zeros((size(FLD.nBot,2)+size(FLD.nTop,2))*1*4*(length(FLD.nBot)-1),4); % [2x2] matrices for 2 surfaces with 1 interacting node
cntApu_O = 0;
a_pu = MP.dZ*FLD.rho*FLD.beta_wBM.*FLD.dcf_32;
% Refer to Pub 1 Supporting Materials.docx Eq. (8)and (9)

switch MP.key_FSI_2Chambers

    case {3}
            nn = size(FLD.nBot,2);
            for in = 1:nn
                % BM
                krow = bcBM;                         % rows
                kcol = (FLD.nBot(:,in)-1)*ndof+fdir;
                [partsApu_O,cntApu_O] = fillT(krow(2:end),kcol(1:end-1),a_pu(2:end)*MP.key_FSI/6/nn,partsApu_O,cntApu_O);
                [partsApu_O,cntApu_O] = fillT(krow(2:end),kcol(2:end),a_pu(2:end)*MP.key_FSI/3/nn,partsApu_O,cntApu_O);
                [partsApu_O,cntApu_O] = fillT(krow(1:end-1),kcol(2:end),a_pu(1:end-1)*MP.key_FSI/6/nn,partsApu_O,cntApu_O);
                [partsApu_O,cntApu_O] = fillT(krow(1:end-1),kcol(1:end-1),a_pu(1:end-1)*MP.key_FSI/3/nn,partsApu_O,cntApu_O);
            end
            
            nn = size(FLD.nTop,2);
            for in = 1:nn
                % RL
                krow = bcTM;                         % rows
                kcol = (FLD.nTop(:,in)-1)*ndof+fdir; % cols
                [partsApu_O,cntApu_O] = fillT(krow(2:end),kcol(1:end-1),-a_pu(2:end)*MP.key_FSI/6/nn,partsApu_O,cntApu_O);
                [partsApu_O,cntApu_O] = fillT(krow(2:end),kcol(2:end),-a_pu(2:end)*MP.key_FSI/3/nn,partsApu_O,cntApu_O);
                [partsApu_O,cntApu_O] = fillT(krow(1:end-1),kcol(2:end),-a_pu(1:end-1)*MP.key_FSI/6/nn,partsApu_O,cntApu_O);
                [partsApu_O,cntApu_O] = fillT(krow(1:end-1),kcol(1:end-1),-a_pu(1:end-1)*MP.key_FSI/3/nn,partsApu_O,cntApu_O);
            end               
end

% Aup = spalloc(Nd.tdof,FLD.tdof,12*FLD.tdof);
if opt_mean_pressure == 0
    surface_interaction = 1;
else
    surface_interaction = 2;
end

partsAup = zeros((size(FLD.nBot,2)+size(FLD.nTop,2))*surface_interaction*4*(length(FLD.nBot)-1),4); % [2x2] matrices for surfaces with 1 interacting node
cntAup = 0;
% tmp = transpose(linspace(0.4,1.2,length(FLD.nBot)));
pAreaB = (MP.dZ*MP.Prop_fit.width_BM(MP.xx(:)).*FLD.dcf_23); % 2 factor is for finding the total width of BM from the half point coordinate
pAreaR = (MP.dZ*MP.Prop_fit.width_BM(MP.xx(:)).*FLD.dcf_23);
% Refer to Pub 1 Supporting Materials.docx Eq.(10) (11) and (12)
switch MP.key_FSI_2Chambers
    
    case {3}
                       
        if opt_mean_pressure == 0
            
            nn = size(FLD.nBot,2);
            for in = 1:nn
                % BM
                krow = (FLD.nBot(:,in)-1)*ndof + fdir;
                kcol1 = bcBM;

                [partsAup,cntAup] = fillT(krow(2:end),kcol1(1:end-1),-pAreaB(2:end)*MP.key_FSI/6/nn,partsAup,cntAup);
                [partsAup,cntAup] = fillT(krow(2:end),kcol1(2:end),-pAreaB(2:end)*MP.key_FSI/3/nn,partsAup,cntAup);
                [partsAup,cntAup] = fillT(krow(1:end-1),kcol1(2:end),-pAreaB(1:end-1)*MP.key_FSI/6/nn,partsAup,cntAup);
                [partsAup,cntAup] = fillT(krow(1:end-1),kcol1(1:end-1),-pAreaB(1:end-1)*MP.key_FSI/3/nn,partsAup,cntAup);

            end
            
            nn = size(FLD.nTop,2);
            for in = 1:nn
                krow = (FLD.nTop(:,in)-1)*ndof + fdir;
                kcol1 = bcTM;

                [partsAup,cntAup] = fillT(krow(2:end),kcol1(1:end-1),pAreaR(2:end)*MP.key_FSI/6/nn,partsAup,cntAup);
                [partsAup,cntAup] = fillT(krow(2:end),kcol1(2:end),pAreaR(2:end)*MP.key_FSI/3/nn,partsAup,cntAup);
                [partsAup,cntAup] = fillT(krow(1:end-1),kcol1(2:end),pAreaR(1:end-1)*MP.key_FSI/6/nn,partsAup,cntAup);
                [partsAup,cntAup] = fillT(krow(1:end-1),kcol1(1:end-1),pAreaR(1:end-1)*MP.key_FSI/3/nn,partsAup,cntAup);

            end

        else
            nn = size(FLD.nBot,2);
            for in = 1:nn
                % BM
                krow = (FLD.nBot(:,in)-1)*ndof + fdir; % column index ; plus one is to shift the node address forward to exclude the first node from FSI BC, becuz it is not a FSI BC in the fluid either
                kcol1 = bcBM;
                [partsAup,cntAup] = fillT(krow(2:end),kcol1(1:end-1),-0.5*pAreaB(2:end)*MP.key_FSI/6/nn,partsAup,cntAup);
                [partsAup,cntAup] = fillT(krow(2:end),kcol1(2:end),-0.5*pAreaB(2:end)*MP.key_FSI/3/nn,partsAup,cntAup);
                [partsAup,cntAup] = fillT(krow(1:end-1),kcol1(2:end),-0.5*pAreaB(1:end-1)*MP.key_FSI/6/nn,partsAup,cntAup);
                [partsAup,cntAup] = fillT(krow(1:end-1),kcol1(1:end-1),-0.5*pAreaB(1:end-1)*MP.key_FSI/3/nn,partsAup,cntAup);

                kcol1 = bcTM;
                [partsAup,cntAup] = fillT(krow(2:end),kcol1(1:end-1),0.5*pAreaB(2:end)*MP.key_FSI/6/nn,partsAup,cntAup);
                [partsAup,cntAup] = fillT(krow(2:end),kcol1(2:end),0.5*pAreaB(2:end)*MP.key_FSI/3/nn,partsAup,cntAup);
                [partsAup,cntAup] = fillT(krow(1:end-1),kcol1(2:end),0.5*pAreaB(1:end-1)*MP.key_FSI/6/nn,partsAup,cntAup);
                [partsAup,cntAup] = fillT(krow(1:end-1),kcol1(1:end-1),0.5*pAreaB(1:end-1)*MP.key_FSI/3/nn,partsAup,cntAup);

            end
            
            nn = size(FLD.nTop,2);
            for in = 1:nn
                krow = (FLD.nTop(:,in)-1)*ndof + fdir; % column index
                kcol1 = bcTM;
                [partsAup,cntAup] = fillT(krow(2:end),kcol1(1:end-1),0.5*pAreaR(2:end)*MP.key_FSI/6/nn,partsAup,cntAup);
                [partsAup,cntAup] = fillT(krow(2:end),kcol1(2:end),0.5*pAreaR(2:end)*MP.key_FSI/3/nn,partsAup,cntAup);
                [partsAup,cntAup] = fillT(krow(1:end-1),kcol1(2:end),0.5*pAreaR(1:end-1)*MP.key_FSI/6/nn,partsAup,cntAup);
                [partsAup,cntAup] = fillT(krow(1:end-1),kcol1(1:end-1),0.5*pAreaR(1:end-1)*MP.key_FSI/3/nn,partsAup,cntAup);

                kcol1 = bcBM;
                [partsAup,cntAup] = fillT(krow(2:end),kcol1(1:end-1),-0.5*pAreaR(2:end)*MP.key_FSI/6/nn,partsAup,cntAup);
                [partsAup,cntAup] = fillT(krow(2:end),kcol1(2:end),-0.5*pAreaR(2:end)*MP.key_FSI/3/nn,partsAup,cntAup);
                [partsAup,cntAup] = fillT(krow(1:end-1),kcol1(2:end),-0.5*pAreaR(1:end-1)*MP.key_FSI/6/nn,partsAup,cntAup);
                [partsAup,cntAup] = fillT(krow(1:end-1),kcol1(1:end-1),-0.5*pAreaR(1:end-1)*MP.key_FSI/3/nn,partsAup,cntAup);
            end
        end
        
    otherwise
        
end

partsAup = partsAup(any(partsAup(:,[3,4]),2),:);
Aup = sparse(partsAup(:,1),partsAup(:,2),partsAup(:,3),Nd.tdof,FLD.tdof);
Apu_O = sparse(partsApu_O(:,1),partsApu_O(:,2),partsApu_O(:,3),FLD.tdof,Nd.tdof);
Aup = -Aup(logical(Nd.mapF2R),FLD.BC);
Apu_O = -Apu_O(FLD.BC,logical(Nd.mapF2R));

AMat.Aup = Aup;
AMat.Apu = Apu_O; % The O comes from the FDomain convention where freq is multiplied later in A1

end