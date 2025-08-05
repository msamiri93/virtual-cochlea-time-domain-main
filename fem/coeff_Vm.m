function [alpha0,Cm0,M] = coeff_Vm(M,Vm0)

%     a = 37.4e-3; % e/kT (1/mV),Note the unit conversion factor 1e-3, V -> mV
%     epsilon = exp(-M.zC*a*(Vm0-M.VC05));
% 
%     Cm0 = M.Qmax*(M.zC*a)*epsilon/(1+epsilon)^2 + M.CLin;
% 
%     alpha0 = -Vm0*M.Qmax*(M.zC*a)^2*(epsilon/(1+epsilon)^2 - 2*epsilon^2/(1+epsilon)^3);
% ****************************************************************************************
% Unit conversion
% Cm: pF
% Gm: nS
% Q: pC

    a = 37.4; % e/kT (1/V),
    epsilon = exp(-M.zC*a*(Vm0-M.VC05)*1e-3);%Note the unit conversion factor 1e-3, V -> mV

    % Capacitance as a function of Vm: Santos-Sacchi, Shen, Zheng, and Dallos 2001
    % 3 times was added to account for 3 OHCs. in the
    % organ_Corti_parameters_wgrad function the 3x was included but later
    % in the code the M structure from OHC_sg was reused that missed the 3x
    % while making the conductance and capacitance values
    Cm0 = (M.Qmax*(M.zC*a)*epsilon/(1+epsilon)^2 + M.CLin);
    % Alpha is Vm*(dCm/dVm) in Yanju Liu Thesis
    alpha0 = (-1e-3*Vm0*M.Qmax*(M.zC*a)^2*(epsilon/(1+epsilon)^2 - 2*epsilon^2/(1+epsilon)^3));%Note the unit conversion factor 1e-3, V -> mV

end