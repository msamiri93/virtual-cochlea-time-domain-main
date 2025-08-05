function [b, pre_in] = evaluate_RHS(time,b0,MP,Nd,El,FLD,OHC,G,R0,R,U,a,v,x,dpo,Vnd)

    opt_cfld = MP.opt_cfld;

    b = b0;

    pre_in = zeros(FLD.tdof,1);

    if MP.vMC == 1
        ind_pre = strncmp(FLD.name,'BCS',3);
    else
        ind_pre = strncmp(FLD.name,'OW',2);
    end

    % initial coding for pure tone only
    ww = 2*pi*MP.freq;
    %%
    if (MP.Estim == 1)
        pre_in(ind_pre,1) = 0;
    else        
        % % half-sin initial
        % raise_time = 5*MP.period;
        % if time < raise_time
        %     h = 0.5*(1 - cos(2*pi*time/(2*raise_time)));
        % % elseif time > MP.T - MP.period
        % %     h = 0.5*(1 + cos(2*pi*time/(2*MP.period)));
        % else
        %     h = 1;
        % end
        pre_ref = 20*1e-6;
        amplitude = pre_ref*10.^(MP.Pstim/20);
        %
        switch MP.stim_type
            case 'random'
                pre_in(ind_pre,1) = ones(sum(ind_pre),1)*stimulationFunction(amplitude, time, freq=MP.freq, type = 'random'); % prescribed pressure applied
            case 'tone'
                pre_in(ind_pre,1) = ones(sum(ind_pre),1)*stimulationFunction(amplitude, time, freq=MP.freq, type = 'tone');
            case 'impulse'
                pre_in(ind_pre,1) = ones(sum(ind_pre),1)*stimulationFunction(amplitude, time, type = 'impulse');
        end     
    end
    RHS = -FLD.L*pre_in;
    row = MP.dof(1).pdof;
    b(row,1) = b(row,1) + RHS;
    pre_in = pre_in(ind_pre,1);
    %%
    % %% F_OHC artificial
    % 
    % % OHC force implementation, proportional to Vm
    % ndof = 6;
    % nOHC = find(strcmpi('OHC',El.name));
    % nOHC = nOHC([OHC.idx]);
    % nd1OHC = El.Nd1(nOHC);
    % nd2OHC = El.Nd2(nOHC);
    % fdir = El.dir(nOHC,:);
    % 
    % fOHC = 0.01*gaussmf(MP.xx,[0.1,2])*sin(ww*time);
    % 
    % for idd = 1:size(fdir,2)
    % 
    %     ff = fdir(:,idd).*(fOHC.');
    % 
    %     rdof = (nd1OHC-1)*ndof + idd;
    %     b(rdof,1) = b(rdof,1) + ff;
    % 
    %     rdof = (nd2OHC-1)*ndof + idd;
    %     b(rdof,1) = b(rdof,1) - ff;
    % 
    % end
    % 
    % %%

    b(MP.BC,1) = b(MP.BC,1) + R0*U(MP.BC,1);    

    row = MP.dof(2).udof;
    b(row(Nd.BC),1) = b(row(Nd.BC),1) + R.Auu_x_a*a(Nd.BC,1) + R.Auu_x_v*v(Nd.BC,1) + R.Auu_x_x*x(Nd.BC,1) - G.Auo*dpo;

    row = MP.dof(4).adof;
    b(row,1) = b(row,1) + R.Aau_x_a*a(Nd.BC,1) + R.Aau_x_v*v(Nd.BC,1) + R.Aau_x_x*x(Nd.BC,1);

    row = MP.dof(3).edof;
    if MP.nonlin
        Aeo = make_Aeo(MP, OHC, Vnd);
        b(row,1) = b(row,1) - Aeo*dpo;
    end

    if opt_cfld == 1
        idx = 1:sum([MP.dof(1:6).n]);
        BC_s = MP.BC(idx,1);
        V = U;
        V(MP.dof(2).udof,1) = v;
        b(MP.dof(7).ldof,1) = b(MP.dof(7).ldof,1) + R.Acm_x_a*U(BC_s,1) + R.Acm_x_v*V(BC_s,1);
    end
end