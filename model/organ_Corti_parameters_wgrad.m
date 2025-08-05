function [StV,SL,OoC] = organ_Corti_parameters_wgrad(MP,nz)
    % d=0.1;    %source: Mistrik 2009
    % BM.x=yBM;    % Basilar membrane displacement (input)
    
    % Stria Vascularis(StV)
    if MP.vMC == 1
        StV_G = 1e4;
    else
        StV_G = 1;
    end

    StV.C = 1/4*(66.7)*ones(nz,1);  % Stria Vascularis capacitance in PF, chosen from Mistrik 2009,modified for 40 micron instead of 60 micron(multiplied by 1.5)
    
    StV.G1 = StV_G*1/4*(6.67*1e3*5)*ones(nz,1);     
    StV.G2 = StV_G*(1.7*1e4)*ones(nz,1); %%%%changed on 08/14/19
        
    StV.E= MP.EP*ones(nz,1);       
    
    % Spiral Limbus (SL)
    SL.G1 = 0.25*(99.5)*ones(nz,1); % Spiral limbus conductance in nS, source : Mistrik 2009,modified for 40 micron instead of 60 micron(multiplied by 1.5)
    SL.G2 = 0.25*(33.33)*ones(nz,1);     % Spiral limbus connectivity conductance in nS, source : Mistrik 2009,modified for 40 micron instead of 60 micron(multiplied by 1.5)
    % organ of Corti (OoC)
    OoC.G1 = 0.25*6.67*1e4*ones(nz,1);   % Organ_of_corti conductance in nS, source: Mistrik 2009,modified for 40 micron instead of 60 micron(multiplied by 1.5) %MSA: temporarily changed conductance by a factor to later figure out why it is so high in the model
    OoC.G2 = 0.25*6.67*1e4*ones(nz,1); % Organ_of_corti connectivity conductance in nS, source: same as OoC.G1 for Mistrik 2009,modified for 40 micron instead of 60 micron(multiplied by 1.5)
end
