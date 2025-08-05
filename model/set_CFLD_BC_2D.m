function CFLD = set_CFLD_BC_2D(CFLD)

    Nd = CFLD.Nd;
    
    Radi = false(Nd.N,1);
    Radi(Nd.ind_Radi) = true;
    % iuFixed = (Radi | Nd.ind_Apex);
    iuFixed = Radi;
    ivFixed = (Nd.ind_Apex | Nd.ind_Base | Nd.ind_Symm) & (~Radi);
    ipFixed = (Nd.ind_Apex | Nd.ind_Base);
%     ipFixed = Nd.ind_Base;
    iwFixed = Nd.ind_Symm;

    CFLD.tdof = Nd.NV*Nd.N;
    a = ones(CFLD.tdof,1);
    tmp = [iuFixed'; ivFixed'; ipFixed'; iwFixed'];
    a(tmp(:)) = 0;         
    CFLD.BC = logical(a);     % boundary condition, indices of fixed boundary DOF == 0
    CFLD.rdof = sum(CFLD.BC);

    Nd.iuFixed = iuFixed;
    Nd.ivFixed = ivFixed;
    Nd.ipFixed = ipFixed;
    Nd.iwFixed = iwFixed;

end