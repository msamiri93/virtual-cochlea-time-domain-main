function AMat = make_AMatrix_FLD_CFLD(MP,FLD,AMat)

    CFLD = FLD.CFLD;
    nd = 2; %number of elements on the edge

    cpdof = ((1:CFLD.Nd.N)-1)*CFLD.Nd.NV + 3; %pressure dof of Corti
    [~,idx] = sort(CFLD.Nd.x(CFLD.Nd.ind_Radi,1));
    cpdof_FSI = cpdof(CFLD.Nd.ind_Radi(idx));
    bcBM = FLD.indBM;
    [~,idx] = sort(FLD.Nd.x(bcBM,1),'ascend');
    bcBM = bcBM(idx);
    nel = FLD.El.N; %number of elements for scalae fluid
    partsMPC = zeros(nd*nd*nel,4);
    cntPC = 0;

    for ie = 1:nel
        knodes = FLD.El.node(ie,:);
        wrap = [knodes,knodes(1)];        
        for id = 1:FLD.El.type               
            j1 = wrap(id); j2 = wrap(id+1);
            i1 = find(bcBM == j1,1);
            i2 = find(bcBM == j2,1);
            if (~isempty(i1) && ~isempty(i2))
                                        
                xe1 = FLD.Nd.x(j1,1); ye1 = FLD.Nd.x(j1,2);
                xe2 = FLD.Nd.x(j2,1); ye2 = FLD.Nd.x(j2,2);
                edge = sqrt((xe2-xe1)^2+(ye2-ye1)^2);
                kk = perm_scale(MP,mean(FLD.CFLD.Nd.Z([i1,i2])));
                width = mean(FLD.CFLD.Nd.width([i1,i2]));
                eqm = -FLD.rho*kk*width*eqm2(edge); % % (-1 normal unit inward)*( + pressure node on bottom of BM)        
                row = [j1,j2]; % should be 1x2
                col = cpdof_FSI([i1,i2]); % should be 1x2
                [partsMPC,cntPC] = fillR(row, col, eqm, partsMPC, cntPC);
            end
        end
    end
    
    partsMPC = partsMPC(any(partsMPC(:,[3,4]),2),:);
    Qpc_O = sparse(partsMPC(:,1), partsMPC(:,2), partsMPC(:,3) + 1i*partsMPC(:,4), FLD.tdof, CFLD.tdof);
    AMat.Qpc_O = Qpc_O(FLD.BC,CFLD.BC);
end

function eqm = eqm2(ds)

    % page 303 Pozrikidis, intro to finite and spectral element methods using matlab
    eqm(1,1) = ds/3;
    eqm(1,2) = ds/6;
    eqm(2,1) = ds/6;
    eqm(2,2) = ds/3;
    
end
