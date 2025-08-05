function S = power_analysis(MP,Nd,El,AMat,FLD,OHC,R,U,uu)

    udof = reshape(1:MP.dof(2).n,6,Nd.Nr,Nd.Nc);

    tt = linspace(0,MP.period,MP.steps+1);
    uf = (2/MP.steps)*uu*exp(-1i*2*pi*MP.freq*tt');

    ww = 2*pi*MP.freq;
    vf = 1i*ww*uf;
    af = ((1i*ww)^2)*uf;

    z0 = zeros(Nd.tdof,1);
    S = struct('kk',z0,'cc',z0,'mm',z0,'ff_met',z0,'ff_ohc',z0,'ff_fld',z0,'ff_cfld',z0);
    % structure
    S.kk(Nd.BC,1) = 0.5*conj(vf(Nd.BC,1)).*AMat.K*uf(Nd.BC,1);
    S.cc(Nd.BC,1) = 0.5*conj(vf(Nd.BC,1)).*AMat.C*vf(Nd.BC,1);
    S.mm(Nd.BC,1) = 0.5*conj(vf(Nd.BC,1)).*AMat.M*af(Nd.BC,1);

    %MET force
    odof = MP.dof(5).odof;
    HB = [OHC.HB]; po0 = [HB.po0];
    po = U(odof,:);
    dpo = po - po0.';
    dpo = (2/MP.steps)*dpo*exp(-1i*2*pi*MP.freq*tt');
    S.ff_met(Nd.BC,1) = 0.5*conj(vf(Nd.BC,1)).*AMat.Auo*dpo;

    %OHC force
    edof = MP.dof(3).edof;
    dve = U(edof,:) - [OHC.Vnd].';
    dve = (2/MP.steps)*dve*exp(-1i*2*pi*MP.freq*tt');
    S.ff_ohc(Nd.BC,1) = 0.5*conj(vf(Nd.BC,1)).*AMat.Aue*dve;

    % fluid force
    pdof = MP.dof(1).pdof;
    ps = U(pdof,:);
    ps = (2/MP.steps)*ps*exp(-1i*2*pi*MP.freq*tt');
    S.ff_fld(Nd.BC,1) = 0.5*conj(vf(Nd.BC,1)).*AMat.Aup*ps(FLD.BC);

    %Corti fluid force
    cdof = MP.dof(6).cdof;
    pc = U(cdof,:);
    pc = (2/MP.steps)*pc*exp(-1i*2*pi*MP.freq*tt');
    S.ff_cfld(Nd.BC,1) = 0.5*conj(vf(Nd.BC,1)).*AMat.Auc*pc(FLD.CFLD.BC);

    fname = fieldnames(S);
    for ii = 1:numel(fname)
        S.(fname{ii}) = sum(reshape(S.(fname{ii}),6,Nd.Nr,Nd.Nc),[1,3]);
    end

    S = struct2table(S);
end