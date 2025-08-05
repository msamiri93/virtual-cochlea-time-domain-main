function M = compute_membrane_capacitance(M, dT)
    
    % MP stands for membrane potential
    % M is a structure that contains the membrane electrical properties
    a = 37.4; % e/kT (1/V)
    
    epsilon = exp(-M.zC*(M.Vm-M.VC05)*a*1e-3);    % Note the unit conversion factor 1e-3, V -> mV
        
    M.Cm = M.Qmax*(M.zC*a)*epsilon/(1+epsilon)^2 + M.CLin; 
    M.Cm = M.Cm; % testing lower capacitance
    M.alpha = 1e-3*M.Vm*M.Qmax*(M.zC*a)^2*(epsilon/(1+epsilon)^2 - 2*epsilon^2/(1+epsilon)^3);
    % Note the unit conversion factor 1e-3, [mV][pC][1/V^2]*1e-3 = [pF]
    
    d_zeta = (dT/M.tauK)*(1/(1+exp(-(M.Vm-M.VG05)/M.eK)) - M.zeta);
    M.zeta = M.zeta + d_zeta;
    
    M.Gm = M.Gmax*M.zeta;    
    
end