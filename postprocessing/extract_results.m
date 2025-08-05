function R = extract_results(MP,Nd,El,U,uu,R)

    nAF = strcmp(Nd.name,'AF');
    nBB = strcmp(Nd.name,'BB');
    nDD = strcmp(Nd.name,'DD');
    nANK = strcmp(El.name,'ANK');
    nd1 = El.Nd1(nANK); nd2 = El.Nd2(nANK);
    odof = MP.dof(5).odof;

    tt = linspace(0,MP.period,MP.steps+1);
    uf = (2/MP.steps)*uu*exp(-1i*2*pi*MP.freq*tt');

    R.uf = reshape(uf,6,Nd.N);

    R.yBM = xx(2,nAF,:);
    R.yRL = xx(2,nBB,:);
    R.xDC = xx(1,nDD,:);
    R.yDC = xx(2,nDD,:);

    R.SPL = MP.Pstim - 25;
    R.sensitive_cochlea = MP.sensitive_cochlea;
    R.freq = MP.freq;
    R.xx = MP.xx;

end