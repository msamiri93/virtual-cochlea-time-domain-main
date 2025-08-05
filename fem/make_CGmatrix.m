function [C,G] = make_CGmatrix(MP,nv,nz,V0,StV,OHC,IHC,SL,OoC)

    partsC = zeros(nv*nv*nz,4);
    cntC = 0;

    partsG = zeros(nv*nv*nz,4);
    cntG = 0;

    for iz=1:nz
        % Capacitance matrix per section
        M = OHC(iz).M;
        Vm0 = V0((iz-1)*5+2)-V0((iz-1)*5+3);
        [a0,Cm0] = coeff_Vm(M,Vm0);

        OHC_Cm = OHC(iz).count*(Cm0 + a0);
        OHC_Cs = OHC(iz).count*OHC(iz).S.C;

        Cm=[OHC_Cs+StV.C(iz)+IHC.Cs(iz), -OHC_Cs, 0, -IHC.Cs(iz), 0;
            -OHC_Cs, OHC_Cs + OHC_Cm, -OHC_Cm, 0, 0;
            0, -OHC_Cm, OHC_Cm, 0, 0;
            -IHC.Cs(iz), 0, 0, IHC.Cs(iz) + IHC.Cm(iz), -IHC.Cm(iz);
            0, 0, 0, -IHC.Cm(iz), IHC.Cm(iz)];

        rdof = 5*iz-4:5*iz;
        [partsC,cntC] = fillR(rdof,rdof,Cm,partsC,cntC);
        
        % Conductance connectivity matrix, off-diagonal, per section
        Gm = [-StV.G2(iz), 0, 0;
                0,-OoC.G2(iz),0;
                0, 0,-SL.G2(iz)];
        i=5*iz-4;
        rdof = [i,i+2,i+4];
        if iz~=1
            cdof = [i-5, i-3, i-1];
            [partsG,cntG] = fillR(rdof,cdof,Gm,partsG,cntG);
        end
        if iz~=nz
        cdof = [i+5, i+7, i+9];
            [partsG,cntG] = fillR(rdof,cdof,Gm,partsG,cntG);
        end

        % Conductance matrix per section
        [Gm0,sGm] = coeff_OM(M,Vm0);
%         az = Gm0 + sGm*75; % MSA: hard-coded value!!!!
        OHC_Gm = OHC(iz).count*Gm0; % + sGm*OHC.Ek(iz); MSA: there is no reference to how the second term is evaluated
        % if ~MP.nonlin
        %     OHC_Gs = OHC(iz).count*OHC(iz).S.Gmax*OHC(iz).ini.po;
        % else
        %     OHC_Gs = 0;
        % end
        OHC_Gs = OHC(iz).count*OHC(iz).S.Gmax*OHC(iz).ini.po;
        StV_t = StV.G1(iz)+StV.G2(iz);
        OoC_t = OoC.G1(iz) + OoC.G2(iz);
        SL_t = SL.G1(iz)+SL.G2(iz);
        G1=[StV_t + OHC_Gs+ IHC.Gs(iz), -OHC_Gs, 0, -IHC.Gs(iz), 0];
        G2=[-OHC_Gs, OHC_Gs + OHC_Gm, -OHC_Gm, 0, 0];
        G3=[0, -OHC_Gm, OoC_t + OHC_Gm, 0, 0];
        G4=[-IHC.Gs(iz), 0, 0, IHC.Gm(iz) + IHC.Gs(iz), -IHC.Gm(iz)];
        G5=[0, 0, 0, -IHC.Gm(iz), SL_t + IHC.Gm(iz)];
        Gm=[G1;G2;G3;G4;G5];

        if ~(iz==1 || iz==nz)
            Gm = Gm + diag([StV.G2(iz) 0 OoC.G2(iz) 0 SL.G2(iz)]);
        end

        rdof = (5*iz-4):5*iz;
        [partsG,cntG] = fillR(rdof,rdof,Gm,partsG,cntG);

    end

    C = sparse(partsC(:,1),partsC(:,2),partsC(:,3),nv*nz,nv*nz);
    G = sparse(partsG(:,1),partsG(:,2),partsG(:,3),nv*nz,nv*nz);
    
end