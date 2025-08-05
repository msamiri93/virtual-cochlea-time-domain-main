function b = initiate_RHS(MP,Nd,OHC,AMat,Ge,V0,Ie)

    HB = [OHC.HB];
    po0 = transpose([HB.po0]);
    V_lin = AMat.Aeo*po0;
    % V0 = Ge\Ie;

    b = zeros(MP.tdof,1);

    row = MP.dof(2).udof;
    b(row(Nd.BC),1) = AMat.Aue*V0;

    row = MP.dof(3).edof;
    
    % b(row,1) = Ie + V_lin;
    if MP.nonlin
        b(row,1) = Ie;
    else
        b(row,1) = Ie + V_lin;
    end

    row = MP.dof(5).odof;
    b(row,1) = po0;
    
end