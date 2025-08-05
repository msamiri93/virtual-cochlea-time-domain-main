function U = initiate_U(MP,OHC)

    HB = [OHC.HB];
    po0 = transpose([HB.po0]);
    Vnd0 = transpose([OHC.Vnd0]);

    U = zeros(MP.tdof,round(MP.period/MP.dt) + 1);

    row = MP.dof(3).edof;
    U(row,1) = Vnd0;

    row = MP.dof(5).odof;
    U(row,1) = po0;

end