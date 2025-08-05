function [OHC,IHC] = model_OHC_2state_circuit(MP, El)

    opt_plot = 0; % 1: plot the convergence

    if isfield(MP,'xx')
        nOHC = length(MP.xx);
    else
        nOHC = length(find(strcmp(El.name,'OHC'))); % number of OHC
    end

    OHC(nOHC) = struct('HB',[],'S',[],'C',[]);

    zz = (-0.5*MP.length_BM:MP.dZ:0.5*MP.length_BM)';
    idx = knnsearch(MP.zz(:),zz(:));
    zz = zz*MP.opt_lgradient;        
    za = 10; zb = 2;        
    if MP.loc<=6
        zz = zz + 1e3*(MP.loc-zb);
    else
        zz = zz + 1e3*(MP.loc-za);
    end

    nv = 5;

    for iz = 1:nOHC

        [S, M] = OHC_electrical_params2(MP,zz(iz));

        HB = HB_transduction_params2(MP, zz(iz));
        HB = ini_OHC(S,HB);
        
        OHC(iz).count = 3;
        OHC(iz).HB = HB;
        OHC(iz).S = S;
        OHC(iz).M = M;
    
        OHC(iz).S.po0 = HB.po0;
        OHC(iz).ini.po = HB.po0;
        OHC(iz).ini.Gs = OHC(iz).ini.po*OHC(iz).S.Gmax;
        OHC(iz).ini.zeta = 0.8;
        OHC(iz).M.zeta = 0;
        OHC(iz).idx = idx(iz);
    end

    IHC = model_IHC(nOHC);

    cnt = 1;
    err = 1;
    if eq(opt_plot,1)
        figure(11);
        h1 = animatedline;
    end

    if MP.sensitive_cochlea
        tol = 1e-12;
    else
        tol = 1e-3;
    end

    while err > tol
        [~,~,~,Vnd0] = mElectrical2_cmod_MSA(MP,OHC,IHC);
        M = [OHC.M];
        ini = [OHC.ini];
        Vnd0 = reshape(Vnd0,nv,nOHC);
        Vm0 = Vnd0(2,:) - Vnd0(3,:);
        zeta0 = (1+exp(-(Vm0-[M.VG05])./[M.eK])).^(-1);
        err = norm(zeta0 - [ini.zeta])/norm([ini.zeta]);

        if eq(opt_plot,1)
            addpoints(h1,cnt,err); drawnow;
        end
       
        for iz = 1:nOHC
            OHC(iz).M.zeta = zeta0(iz);
            OHC(iz).ini.zeta = zeta0(iz);
        end

        cnt = cnt + 1;
        if cnt == 500
            num_s = num2str(cnt,3);
            err_s = ['Solution did nor converge in ',num_s,' iterations.'];
            error('myApp:argChk',err_s)
        end
    end

    for iz = 1:nOHC
        M = OHC(iz).M;
        
        % Initial value of OHC membrane potential is given assuming
        % 0.45 resting open probability
        OHC(iz).Vnd0 = Vnd0((iz-1)*nv+1:iz*nv);
        OHC(iz).Vnd = OHC(iz).Vnd0;
        M.Vm0 = Vnd0((iz-1)*nv+2) - Vnd0((iz-1)*nv+3);
        M.Vm = M.Vm0;
        M = compute_membrane_capacitance(M, 0);

        OHC(iz).M = M;
        OHC(iz).ini.Gm = M.Gm;
        OHC(iz).ini.Cm = M.Cm;
        OHC(iz).ini.alpha = M.alpha;
        OHC(iz).ini.zeta = M.zeta;
    end

end