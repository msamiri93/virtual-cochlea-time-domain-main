function FLD = make_CFLD_matrix(MP,FLD)

    CFLD = FLD.CFLD;
    MP.rho = CFLD.rho;
    El = CFLD.El;
    Nd = CFLD.Nd;

    NV = Nd.NV; %number of variables
    N = El.type*NV; % number of dof for an element

    partsMA = zeros(N*N*El.N,4); cntA = 0;
    partsMR = zeros(N*N*El.N,4); cntR = 0;
    b = zeros(NV*Nd.N,1); % right hand side

    for ie = 1:El.N
        knodes = El.node(ie,:);

        [A_e,R_e,F_e] = element_matrix(ie,El,Nd,MP);
        [partsMA, cntA] = fillM(NV, knodes', A_e, partsMA, cntA);
        [partsMR, cntR] = fillM(NV, knodes', R_e, partsMR, cntR);
        b = fillV(NV, knodes', F_e, b);
                         
    end

    Acc = sparse(partsMA(:,1), partsMA(:,2), partsMA(:,3), CFLD.tdof, CFLD.tdof);
    Rcc = sparse(partsMR(:,1), partsMR(:,2), partsMR(:,3), CFLD.tdof, CFLD.tdof);

    cdof = CFLD.BC;

    CFLD.Acc = Acc(cdof,cdof);
    CFLD.Rcc = Rcc(cdof,cdof);
    CFLD.b = b;
    FLD.CFLD = CFLD;
end

function [A_e,R_e,F_e] = element_matrix(ie,El,Nd,MP)

    [~, ~, w] = set_Gauss_local_variables(2,2,El.NQ);
    NQ = El.NQ;
    dt_inv = 1/MP.dt;
    nu = MP.nu;
    NV = Nd.NV;
    rho_inv = 1/MP.rho;
    A_e = zeros(NV*El.type);
    R_e = zeros(NV*El.type);
    F_e = zeros(NV*El.type,1);
    L = zeros(4,NV*El.type); %first dim is number of equations
    R = zeros(4,NV*El.type); %first dim is number of equations
    hs = El.elm_hs(:,:,ie);


    for iq = 1:NQ
        
        cf = hs(iq)*w(iq);

        A1 = [0 0 (1-MP.af)*rho_inv 0;0 0 0 -(1-MP.af)*nu; 0 -1 0 0; 1 0 0 0];
        A2 = [0 0 0 (1-MP.af)*nu;0 0 (1-MP.af)*rho_inv 0; 1 0 0 0; 0 1 0 0];
        A0 = [dt_inv 0 0 0; 0 dt_inv 0 0; 0 0 0 1; 0 0 0 0];

        R1 = [0 0 -MP.af*rho_inv 0;0 0 0 MP.af*nu; 0 0 0 0; 0 0 0 0];
        R2 = [0 0 0 -MP.af*nu;0 0 -MP.af*rho_inv 0; 0 0 0 0; 0 0 0 0];
        R0 = [dt_inv 0 0 0; 0 dt_inv 0 0; 0 0 0 0; 0 0 0 0];        

        f = [0; 0; 0; 0];
       
        for in = 1:El.type

            psi = El.elm_psi(in,iq,ie);
            gpsix = El.elm_gpsi(in,iq,ie);
            gpsiy = El.elm_gpsi(in + El.type,iq,ie);

            L(:,NV*(in-1)+1:NV*in) = gpsix*A1 + gpsiy*A2 + psi*A0;
            R(:,NV*(in-1)+1:NV*in) = gpsix*R1 + gpsiy*R2 + psi*R0;
        end
        
        A_e = A_e + (L.')*L*cf;
        R_e = R_e + (L.')*R*cf;
        F_e = F_e + (L.')*f*cf;
        
    end
end