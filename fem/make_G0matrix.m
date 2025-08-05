function [G,f] = make_G0matrix(nv, nz,StV,OHC,IHC,SL,OoC)
    % % This is the analysis of whole three scalae circuit.This is giving output for dynamic condition.
    % % This circuit is for creating the differential equation for the circuit
    % % First written in May,2014 as a part of the course project BME 404
    % % The circuit is a combination Mistrik et al. 2009 paper and NAMLAB circuit
    % % The circuit takes time, voltage and basilar membrane displacement as input and returns differential voltage as output

    partsG = zeros(nv*nv*nz,4);
    cntG = 0;
    f = zeros(nv*nz,1);
 
    for iz=1:nz
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
        OHC_Gs = OHC(iz).count*OHC(iz).S.Gmax*OHC(iz).ini.po;
        StV_t = StV.G1(iz)+StV.G2(iz);
        OHC_Gm = OHC(iz).count*OHC(iz).M.Gmax*OHC(iz).ini.zeta;
        OoC_t = OoC.G1(iz) + OoC.G2(iz);
        SL_t = SL.G1(iz)+SL.G2(iz);

        G1=[StV_t + OHC_Gs + IHC.Gs(iz), -OHC_Gs, 0, -IHC.Gs(iz), 0];
        G2=[-OHC_Gs, OHC_Gs + OHC_Gm, -OHC_Gm, 0, 0];
        G3=[0, -OHC_Gm, OoC_t + OHC_Gm, 0, 0];
        G4=[-IHC.Gs(iz), 0, 0, IHC.Gs(iz) + IHC.Gm(iz), -IHC.Gm(iz)];
        G5=[0, 0, 0, -IHC.Gm(iz), SL_t + IHC.Gm(iz)];
        Gm=[G1;G2;G3;G4;G5];

        if ~(iz==1 || iz==nz)
            Gm = Gm + diag([StV.G2(iz) 0 OoC.G2(iz) 0 SL.G2(iz)]);
        end

        rdof = (5*iz-4):5*iz;
        [partsG,cntG] = fillR(rdof,rdof,Gm,partsG,cntG);

        fs=[StV.E(iz)*StV.G1(iz);-OHC(iz).M.Ek*OHC_Gm;OHC(iz).M.Ek*OHC_Gm;-IHC.Ek(iz)*IHC.Gm(iz);IHC.Ek(iz)*IHC.Gm(iz)];
        f(5*iz-4:5*iz,1)=fs;

    end
    
    G = sparse(partsG(:,1),partsG(:,2),partsG(:,3),nv*nz,nv*nz);    

end