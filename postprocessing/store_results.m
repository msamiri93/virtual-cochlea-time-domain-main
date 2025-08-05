function R = store_results(tag,ct,PoI,MP,fMP,Nd,El,OHC,FLD,AMat,Uf)

    R.freq = fMP.ffsym(fMP.ind);

    R.Uf = Uf(:,fMP.ind);

    % scala fluid
    pdof = MP.dof(1).pdof;
    R.pre = Uf(pdof,fMP.ind);
    
    % structure
    ndof = 6;
    um2nm = 1e3;
    udof = MP.dof(2).udof;
    uu_f = Uf(udof,fMP.ind);
    idx = reshape(transpose((1:Nd.tdof)),ndof,Nd.N);
    dx_f = um2nm*uu_f(idx(1,:),:);
    dy_f = um2nm*uu_f(idx(2,:),:);
    dz_f = um2nm*uu_f(idx(3,:),:);

    R.yBM = dy_f(PoI.zBM,:);
    R.yRL = dy_f(PoI.nRL,:);
    R.xRL = dx_f(PoI.nRL,:);
    R.yTM = dy_f(PoI.zTM,:);
    R.xTM = dx_f(PoI.zTM,:);
    R.aDC = dy_f(PoI.aDC,:);
    
    nANK = find(strcmpi('ANK',El.name));
    nd1ANK = El.Nd1(nANK); 
    nd2ANK = El.Nd2(nANK);
    R.xHB = dx_f(nd2ANK,:) - dx_f(nd1ANK,:);
    
    if MP.key_fOHC || MP.key_fMET % active cochlea condition

        % voltage
        nv = 5;
        edof = MP.dof(3).edof;
        v = inverse_fourier_transform(Uf(edof,:));    
        idx = reshape(transpose(1:MP.dof(3).n),nv,MP.nz);
        R.Vm = v(idx(2,:),:) - v(idx(3,:),:);
        V0 = reshape([OHC.Vnd0],[],1);

        nOHC = find(strcmpi('OHC',El.name));
        fdir = El.dir(nOHC,:);
%         nd1OHC = El.Nd1(nOHC);
        nd2OHC = El.Nd2(nOHC);
        f = zeros(MP.dof(2).n,fMP.L); % full vector
        f(Nd.BC,:) = -AMat.Aue*(v-V0); % update using reduced matrix
        dir = 2;
%         R.fOHC = f(idx(dir,nd2OHC),:)./fdir(:,dir); % you can project back using any of the three components
    end

end