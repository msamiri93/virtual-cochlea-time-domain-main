function [Gm0,sGm] = coeff_OM(M,Vm0)

    zeta0 = (1+exp(-(Vm0-M.VG05)/M.eK))^(-1);
    % 3 times was added to account for 3 OHCs. in the
    % organ_Corti_parameters_wgrad function the 3x was included but later
    % in the code the M structure from OHC_sg was reused that missed the 3x
    % while making the conductance and capacitance values
    Gm0 = (M.Gmax*zeta0);
    sGm = (-M.Gmax*zeta0^2*exp(-(Vm0-M.VG05)/M.eK)*(-1/M.eK));

end