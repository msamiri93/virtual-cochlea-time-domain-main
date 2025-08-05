function HB = ini_OHC(S,HB)
    HB.xs0 = 0;
    HB.xa0 = 0;
    HB.xa = HB.xa0;
    HB.xx = HB.xs0;
    
    dE = HB.z*(HB.xs0 - HB.xa0 - HB.Xo);
    kBT = 4e-3;
    kCO = HB.kF*exp(0.5*dE/kBT);      % rate from closed to open state
    kOC = HB.kR*exp(-0.5*dE/kBT);     % rate from open to closed state
    po0 = kCO/(kCO + kOC);
    HB.po0 = po0;
    HB.po = HB.po0;
    S.Gs0 = S.Gmax*po0;
    HB.xa = HB.xa0;
    HB.xx = HB.xs0;

%     a = 37.4;
%     epsilon = exp(-M.zC*(M.Vm0-M.VG05)*a*1e-3);    % Note the unit conversion factor 1e-3, mV -> V
%     M.Cm0 = M.Qmax*(M.zC*a)*epsilon/(1+epsilon)^2 + M.CLin; 
%     M.alpha0 = 1e-3*M.Vm0*M.Qmax*(M.zC*a)^2*(epsilon/(1+epsilon)^2 - 2*epsilon^2/(1+epsilon)^3);
% 
%     M.Vm = M.Vm0;
%     M.Cm = M.Cm0;
%     M.Gm = M.Gm0;
%     M.zeta = M.zeta0;
%     M.alpha = M.alpha0;


end