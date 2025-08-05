function U = update_U(time,MP,Nd,FLD,u,U,pre_in)

    if MP.vMC == 1
        ind_pre = strncmp(FLD.name,'BCS',3);
    else
        ind_pre = strncmp(FLD.name,'OW',2);
    end

    % initial coding for pure tone only
    ww = 2*pi*MP.freq;
    %
    if (MP.Estim == 1)
        U(MP.dof(1).pdof(ind_pre),1) = 0;
    else
        U(MP.dof(1).pdof(ind_pre),1) = pre_in; % prescribed pressure applied
    end
    
    U(MP.BC) = u;

    sdof = (Nd.Master_Slave(:,2)-1)*MP.NDOF + (1:3);
    mdof = (Nd.Master_Slave(:,1)-1)*MP.NDOF + (1:3);
    udof = MP.dof(2).udof;
    U(udof(sdof),1) = U(udof(mdof),1);

end