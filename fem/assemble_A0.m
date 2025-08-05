function [A0,R0,L,R,G,El,FLD,OHC,Ce,Ge,Ie,V0] = assemble_A0(Nd,El,OHC,IHC,FLD,MP)

    if MP.Visc == 0
        [A0,R0,L,R,G,El,FLD,OHC,Ce,Ge,Ie,V0] = assemble_inviscid(Nd,El,OHC,IHC,FLD,MP);
    end

end

% #########################################################################
% ###################### INVISCID MATRICES ASSEMBLY #######################
% #########################################################################

function [A0,R0,L,R,G,El,FLD,OHC,Ce,Ge,Ie,V0] = assemble_inviscid(Nd,El,OHC,IHC,FLD,MP)
% %
% % Assembe frequency-independent submatrices that do not need to update with stimulating frequencies
% %

    opt_cfld = MP.opt_cfld;

    %reduced degrees of freedom
    pdof = FLD.rdof;
    udof = Nd.rdof;
    edof = MP.dof(3).n;
    adof = MP.dof(4).n;
    odof = MP.dof(5).n;

    G.Apu = sparse(pdof,udof);
    G.Ape = sparse(pdof,edof);
    G.Apa = sparse(pdof,adof);
    G.Apo = sparse(pdof,odof);

    G.Auu = sparse(udof,udof);
    G.Aua = sparse(udof,adof);
    G.Auo = sparse(udof,odof);
    G.Aue = sparse(udof,edof);

    G.Aep = sparse(edof,pdof);
    G.Aeu = sparse(edof,udof);
    G.Aee = sparse(edof,edof);
    G.Aea = sparse(edof,adof);
    G.Aeo = sparse(edof,odof);

    G.Aap = sparse(adof,pdof);
    G.Aae = sparse(adof,edof);
    G.Aau = sparse(adof,udof);
    G.Aaa = sparse(adof,adof);
    G.Aao = sparse(adof,odof);

    G.Aop = sparse(odof,pdof);
    G.Aou = sparse(odof,udof);
    G.Aoe = sparse(odof,edof);
    G.Aoa = sparse(odof,adof);
    G.Aoo = sparse(odof,odof); 
    
    if opt_cfld == 1

        cdof = FLD.CFLD.rdof;

        G.Qpc = sparse(pdof,cdof);
        G.Auc = sparse(udof,cdof); % place holder, should be defined
        G.Acc = sparse(cdof,cdof); % place holder, should be defined

        G.Acu = transpose(G.Auc);        
        G.Aec = sparse(edof,cdof);
        G.Aac = sparse(adof,cdof);
        G.Aoc = sparse(odof,cdof);
    end

    El.mat(El.N) = struct('m',[],'k',[],'c',[]);
    [M,El] = createM(Nd,El);
    [K,El] = createK(Nd,El);
    [C,El] = createC2(Nd,El,MP);
    
    G.M = M;
    G.C = C;
    G.K = K;

    L.Auu = (1-MP.am)*M + (1-MP.af)*(MP.dt*MP.gamma*C + MP.dt^2*MP.beta*K);
    R.Auu_x_a = -( MP.am*M + (1-MP.af)*(MP.dt*(1-MP.gamma)*C + MP.dt^2*(0.5-MP.beta)*K) );
    R.Auu_x_v = - ( C + (1-MP.af)*MP.dt*K );
    R.Auu_x_x = - K; 

    [G, FLD] = make_FE_FLDmatrix(MP,G,FLD,opt_cfld);
    
    L.App = (1-MP.am)*G.App; % App: gpsi*gpsi
    R.App = - MP.am*G.App;

    [G,FLD] = make_AMatrix_FSI_2(G,Nd,El,FLD,MP);

    L.Aup = (1-MP.am)*G.Aup;
    R.Aup = - MP.am*G.Aup;

    L.Apu = (1-MP.am)*G.Apu;
    R.Apu = - MP.am*G.Apu;

    if opt_cfld == 1

        FLD = make_CFLD_matrix(MP,FLD);
        [L.Acm,R.Acm_x_a,R.Acm_x_v] = make_peristalsis_effect(MP, FLD);
        G = make_AMatrix_FSI_CFLD(Nd,FLD,G,MP);
        G = make_AMatrix_FLD_CFLD(MP,FLD,G); % permeability relation between the two fluid spaces

        L.Acc = FLD.CFLD.Acc;
        R.Acc = FLD.CFLD.Rcc;

        L.Auc = (1-MP.am)*G.Auc;
        R.Auc = - MP.am*G.Auc;
    end

    % C and G matrices at equilibrium state
    [Ce,Ge,Ie,V0] = mElectrical2_cmod_MSA(MP,OHC,IHC); 

    L.Aee = Ce/MP.dt + (1-MP.af)*Ge;
    R.Aee = Ce/MP.dt - MP.af*Ge;

    if ~MP.nonlin
        G.Aeo = make_lin_Aeo(MP, OHC);
    end

    L.Aeo = (1-MP.af)*G.Aeo;
    R.Aeo = - MP.af*G.Aeo;
    
    HB = [OHC.HB];
    G.Aau = make_Aau(MP,Nd,El,OHC);
    [G.Aaa,G.Iaa] = make_Aaa(MP,HB);       

    L.Aau = MP.dt^2*MP.beta*G.Aau;
    R.Aau_x_a = - MP.dt^2*(0.5-MP.beta)*G.Aau;
    R.Aau_x_v = - MP.dt*G.Aau;
    R.Aau_x_x = - G.Aau;

    % % rate of hair bundle deflection
    % L.Aau = MP.dt*MP.gamma*G.Aau;    
    % R.Aau_x_a = - MP.dt*(1-MP.gamma)*G.Aau;
    % R.Aau_x_v = - G.Aau;
    % R.Aau_x_x = sparse(size(G.Aau,1),size(G.Aau,2));

    L.Aaa = G.Iaa/MP.dt + G.Aaa;
    R.Aaa = G.Iaa/MP.dt;

    if MP.sensitive_cochlea        
        [G.Auo, ~] = make_Auo(MP,Nd,El,OHC,HB);
        [G.Aue, ~] = make_Aue(MP,Nd,El,OHC);
    end

    % MET force is fully rhs
    alpha = 0;
    L.Auo = 0*((1-MP.af)*G.Auo + alpha*G.Auo/MP.dt);
    R.Auo = 0*(- MP.af*G.Auo + alpha*G.Auo/MP.dt);

    L.Aue = (1-MP.af)*G.Aue;
    R.Aue = - MP.af*G.Aue;
    
    % active comes later, because G.Auo and G.Aue are needed for rhs vector
    % for the zero frequency where po0 and Vm0 are non-zero

    if ~MP.nonlin
    % po wrt xHB and xa: po - d(po)dxHB*xHB - d(po)dxa*xa = po_rest - d(po)dxHB*xHB_rest - d(po)dxa*xa_rest
        G = make_lin_Aou_Aoa(G, MP, Nd, El, OHC); % matrices for d(po)dxHB and d(po)dxa
    end

    L.Aou = (1-MP.af)*G.Aou;
    R.Aou = - MP.af*G.Aou;
    
    L.Aoa = (1-MP.af)*G.Aoa;
    R.Aoa = - MP.af*G.Aoa;

    L.Aoo = (1-MP.af)*G.Aoo;
    R.Aoo = - MP.af*G.Aoo;
    
    if opt_cfld == 1
        A0 = [ L.App, L.Apu, G.Ape, G.Apa, G.Apo, G.Qpc;
               L.Aup, L.Auu, L.Aue, G.Aua, L.Auo, L.Auc;    
               G.Aep, G.Aeu, L.Aee, G.Aea, L.Aeo, G.Aec;
               G.Aap, L.Aau, G.Aae, L.Aaa, G.Aao, G.Aac;
               G.Aop, L.Aou, G.Aoe, L.Aoa, L.Aoo, G.Aoc;
              G.Qpc', G.Acu,G.Aec',G.Aac',G.Aoc', L.Acc];

        R0 = [ R.App, R.Apu, G.Ape, G.Apa, G.Apo, G.Qpc;
               R.Aup, G.Auu, R.Aue, G.Aua, R.Auo, R.Auc;    
               G.Aep, G.Aeu, R.Aee, G.Aea, R.Aeo, G.Aec;
               G.Aap, G.Aau, G.Aae, R.Aaa, G.Aao, G.Aac;
               G.Aop, R.Aou, G.Aoe, R.Aoa, R.Aoo, G.Aoc;
              G.Qpc', G.Acu,G.Aec',G.Aac',G.Aoc', R.Acc];  

        if MP.nonlin  
            odof_r = pdof + udof + edof + adof + (1:odof);
            idx = true(size(A0,1),1);
            idx(odof_r) = false;
            A0 = A0(idx,idx);
            R0 = R0(idx,idx);
        end
        
        ldof = length(MP.dof(7).dof_r);
        dof = size(L.Acm,2);
        LA = sparse(ldof,dof);
        LL = sparse(ldof,ldof);
        A0 = [ A0, L.Acm.';
               L.Acm, LL  ];

        R0 = [ R0, LA.';
           LA, LL  ];
    else
        % A0 = [ G.App, G.Apu, G.Ape, G.Apa, G.Apo;
        %        G.Aup, G.Auu, G.Aue, G.Aua, G.Auo;    
        %        G.Aep, G.Aeu, G.Aee, G.Aea, G.Aeo;
        %        G.Aap, G.Aau, G.Aae, G.Aaa, G.Aao;
        %        G.Aop, G.Aou, G.Aoe, G.Aoa, G.Aoo];
        A0 = [ L.App, L.Apu, G.Ape, G.Apa, G.Apo;
               L.Aup, L.Auu, L.Aue, G.Aua, L.Auo;    
               G.Aep, G.Aeu, L.Aee, G.Aea, L.Aeo;
               G.Aap, L.Aau, G.Aae, L.Aaa, G.Aao;
               G.Aop, L.Aou, G.Aoe, L.Aoa, L.Aoo];

        R0 = [ R.App, R.Apu, G.Ape, G.Apa, G.Apo;
               R.Aup, G.Auu, R.Aue, G.Aua, R.Auo;    
               G.Aep, G.Aeu, R.Aee, G.Aea, R.Aeo;
               G.Aap, G.Aau, G.Aae, R.Aaa, G.Aao;
               G.Aop, R.Aou, G.Aoe, R.Aoa, R.Aoo];

        if MP.nonlin  
            odof_r = pdof + udof + edof + adof + (1:odof);
            idx = true(size(A0,1),1);
            idx(odof_r) = false;
            A0 = A0(idx,idx);
            R0 = R0(idx,idx);
        end
    end

end

