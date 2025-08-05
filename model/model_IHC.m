function IHC = model_IHC(n)
    % % Inner Hair Cell (IHC)
% %     IHC.Cs = ones(n,1);% OHC steriocilial capacitance in pF, source: Mistrik 2009
% %     IHC.Cm = 10*ones(n,1); % OHC basolater linear capacitance in pF, source: Mistrik 2009
% %     IHC.Gs_max = 0.25*(99.5)*ones(n,1); %IHC steriocilial conductance(maximum)in nS, source: Mistrik 2009
% %     IHC.Gm = 0.25*(33.333)*ones(n,1);        %IHC basolateral conductance in nS, source: Mistrik 2009
% %     IHC.Ek = 100*ones(n,1); %IHC equilibrium potential in mV, source: Mistrik 2009
% %     poi = 0.1;
% %     IHC.Gs = IHC.Gs_max*poi;
%     IHC.Cs = 1e-24*ones(n,1);% OHC steriocilial capacitance in pF, source: Mistrik 2009
%     IHC.Cm = 1e24*10*ones(n,1); % OHC basolater linear capacitance in pF, source: Mistrik 2009
%     IHC.Gs_max = 1e-24*0.25*(99.5)*ones(n,1); %IHC steriocilial conductance(maximum)in nS, source: Mistrik 2009, "The deviant factor is 0.6 "
%     IHC.Gm = 1e24*0.25*(33.333)*ones(n,1);        %IHC basolateral conductance in nS, source: Mistrik 2009
%     IHC.Ek = 0*75*ones(n,1);
%     poi = 0.1;
%     IHC.Gs = IHC.Gs_max*poi;

%     IHC.Cs = ones(n,1);% OHC steriocilial capacitance in pF, source: Mistrik 2009
    IHC.Cs = 1e-3*ones(n,1);% OHC steriocilial capacitance in pF, source: Mistrik 2009
%     IHC.Cm = 10*ones(n,1); % OHC basolater linear capacitance in pF, source: Mistrik 2009
    IHC.Cm = 1e-3*10*ones(n,1); % OHC basolater linear capacitance in pF, source: Mistrik 2009
%     IHC.Gs_max = 0.25*(99.5)*ones(n,1); %IHC steriocilial conductance(maximum)in nS, source: Mistrik 2009

    IHC.Gs_max = 1e-3*0.25*(99.5)*ones(n,1); %IHC steriocilial conductance(maximum)in nS, source: Mistrik 2009

%     IHC.Gm = 0.25*(33.333)*ones(n,1);        %IHC basolateral conductance in nS, source: Mistrik 2009

    IHC.Gm = 1e-3*0.25*(33.333)*ones(n,1);        %IHC basolateral conductance in nS, source: Mistrik 2009
    IHC.Ek = 75*ones(n,1); %IHC equilibrium potential in mV, source: Mistrik 2009
    poi = 0.1;
    IHC.Gs = IHC.Gs_max*poi;

end