function  [Nd,El,OHC,IHC,MP,FLD] = initiate(MP)

    global coef
    opt_fin = 0;  % 0: generate new fluid-domain mesh, 1: use existing mesh grid
    opt_MP_set = 1;         % 0: the MP set before 2021.04
    
    dim = 3;    % dimension 2 or 3-D

    opt_lgradient = 1; % longitudinal gradient
    if MP.vMC == 1, opt_lgradient = 1; end % to control the gradient for testing purposes

    load('./assets/geom/GP_rev0425c.mat','GP');
    % load('./geom/GP_tmp.mat','GP');
    Prop = GP;

    load('./assets/geom/2021_model_geometrical_properties_rev0617a.mat','GP');
    % for i = 1:numel(GP) % for now have everything on log space
    %     GP(i).fit_type = 'exp';
    %     if strcmp(GP(i).name,'thick_BMp')
    %         GP(i).fit_type = 'poly3';
    %     end
    % end    
    GP = rmfield(GP,'descripton');
    Prop = [Prop,GP];

    load('./assets/geom/2021_model_elastic_properties.mat','YP');
    for i = 1:numel(YP) % for now have everything on log space
        YP(i).fit_type = 'exp';        
    end
    YP = rmfield(YP,'descripton');
    Prop = [Prop,YP];
    % Prop = apply_model_coefs(Prop);
    Prop_fit = interpolate_properties(Prop,plot=false);
    Prop = apply_model_coefs(Prop);
    coef_fit = interpolate_properties(Prop,plot=false);
    fields = fieldnames(Prop_fit);
    for i = 1:numel(fields)
        Prop_fit.(fields{i}) = @(x) Prop_fit.(fields{i})(x).*coef_fit.(fields{i})(x);
    end
    Prop_fit = eval_damping_coef(Prop_fit);
    MP.Prop_fit = Prop_fit;
    [Nd, El, MP] = model_3D(MP, dim, opt_lgradient, opt_MP_set);
    MP.opt_lgradient = opt_lgradient;

    if MP.vMC == 1
        MP.EP = 0; % endocochlear potential in mV
    else
        MP.EP = 0.5*90;
    end
    
    MP.EK = 75; % cell equilibrium potential in mV,note the sign, MSA: is this even used?!
    Nd = set_boundary_conditions(MP, Nd, opt_MP_set);
    [OHC, IHC] = model_OHC_2state_circuit(MP,El);
    FLD = model_FLD_cochlea(MP,MP.H,Nd);
    if opt_fin == 1
        input = load('./assets/mesh_12mm_cochlea_quad_3.mat','FLD');
        if isfield(input.FLD,'p')
            FLD.h0 = input.FLD.h0;
            FLD.HOC = input.FLD.HOC;
            FLD.L = input.FLD.L;
            FLD.Nd.x = input.FLD.p;
            FLD.El.node = input.FLD.t;
            FLD.El.N = size(FLD.El.node,1);
            FLD.Nd.N = size(FLD.Nd.x,1);
            FLD.El.type = 3;
        end                
    else          
        FLD = FLD_mesh(MP,MP.H,FLD);
        % the following lines for when quadrilateral elements are used
        if eq(FLD.El.type,4) 
            FLD.El.NQ = 4;
%             set_GaussQuad_as_a_global_variable(FLD.El.NQ);
            [xi, eta, ~] = set_Gauss_local_variables(2,2,FLD.El.NQ);
            FLD.El = evaluate_interp_func(FLD.Nd,FLD.El,xi,eta);
        end
    end
    
    if MP.vMC == 1
        FLD = set_FEfluid_BCs_MC(MP,FLD);
    else
        FLD = set_FEfluid_BCs(MP,FLD);
    end

    if MP.opt_cfld == 1
        CFLD = model_CFLD(MP,Nd);
        CFLD = CFLD_mesh_2D(MP,Nd,CFLD);
        
        CFLD.El.NQ = 4;
        CFLD = set_CFLD_BC_2D(CFLD);
        [xi, eta, ~] = set_Gauss_local_variables(2,2,CFLD.El.NQ);
        CFLD.El = evaluate_interp_func(CFLD.Nd,CFLD.El,xi,eta);

        FLD.CFLD = CFLD;
    end

    if MP.opt_cfld == 1      
        % permeability condition
        alpha = MP.p_alpha;
        beta = MP.p_beta;
        ww = 2*pi*MP.freq;
        kk = alpha;% + 1i*beta/ww;
        MP.kk = kk;
    end

    % initiating degree of freedom structure
    
    MP.dof(1).name = 'scala pressure';
    MP.dof(1).n = FLD.tdof;
    MP.dof(1).pdof = 1:MP.dof(1).n;
    MP.dof(1).dof_r = 1:FLD.rdof;
    last_dof_r = MP.dof(1).dof_r(end);

    MP.dof(2).name = 'structural motion';
    MP.dof(2).n = Nd.tdof;
    MP.dof(2).udof = MP.dof(1).n + (1:MP.dof(2).n);
    MP.dof(2).dof_r = last_dof_r + (1:Nd.rdof);
    last_dof_r = MP.dof(2).dof_r(end);

    nz = length(MP.xx);
    if MP.ind_MET == 2
        nv = 5; % 5 for node voltage
        MP.dof(3).name = '5 node circuit';
        MP.dof(3).n = nv*nz;
        MP.dof(3).edof = sum([MP.dof([1,2]).n]) + (1:MP.dof(3).n);
        MP.dof(3).dof_r = last_dof_r + (1:MP.dof(3).n);
        last_dof_r = MP.dof(3).dof_r(end);

        MP.dof(4).name = 'adaptation x_a';
        MP.dof(4).n = nz;
        MP.dof(4).adof = sum([MP.dof((1:3)).n]) + (1:MP.dof(4).n);
        MP.dof(4).dof_r = last_dof_r + (1:MP.dof(4).n);
        last_dof_r = MP.dof(4).dof_r(end);

        MP.dof(5).name = 'channel open probability';
        MP.dof(5).n = nz;
        MP.dof(5).odof = sum([MP.dof((1:4)).n]) + (1:MP.dof(5).n);
        MP.dof(5).dof_r = last_dof_r + (1:MP.dof(5).n);
        last_dof_r = MP.dof(5).dof_r(end);
    else
        error('10 state channel is not implemented in this code');
    end
    
    if MP.opt_cfld == 1
        MP.dof(6).name = 'Corti fluid';
        MP.dof(6).n = CFLD.tdof;
        MP.dof(6).cdof = sum([MP.dof((1:5)).n]) + (1:MP.dof(6).n);
        MP.dof(6).dof_r = last_dof_r + (1:CFLD.rdof);
        last_dof_r = MP.dof(6).dof_r(end);

        MP.dof(7).name = "Peristalsis";
        MP.dof(7).n = CFLD.Nd.nz;
        MP.dof(7).ldof = sum([MP.dof((1:6)).n]) + (1:MP.dof(7).n);
        BC = reshape(FLD.CFLD.BC,4,FLD.CFLD.Nd.N);
        p_BC = BC(2,FLD.CFLD.Nd.ind_Radi);
        ldof_r = sum(BC(2,FLD.CFLD.Nd.ind_Radi));
        MP.dof(7).dof_r = last_dof_r + (1:ldof_r);
    else
        MP.dof(6).n = 0;
        MP.dof(7).n = 0;
    end

    if MP.opt_cfld == 1
        if ~MP.nonlin
            MP.BC = [FLD.BC; Nd.BC; true(MP.dof(3).n,1); true(MP.dof(4).n,1);...
                    true(MP.dof(5).n,1); CFLD.BC; p_BC(:)];
        else
            MP.BC = [FLD.BC; Nd.BC; true(MP.dof(3).n,1); true(MP.dof(4).n,1);...
                    false(MP.dof(5).n,1); CFLD.BC; p_BC(:)];
        end
        MP.BC_psv = [FLD.BC; Nd.BC; false(MP.dof(3).n,1); false(MP.dof(4).n,1);...
        false(MP.dof(5).n,1); CFLD.BC; p_BC(:)];
    else
        if ~MP.nonlin
            MP.BC = [FLD.BC; Nd.BC; true(MP.dof(3).n,1); true(MP.dof(4).n,1);...
                    true(MP.dof(5).n,1)];
        else
            MP.BC = [FLD.BC; Nd.BC; true(MP.dof(3).n,1); true(MP.dof(4).n,1);...
                    false(MP.dof(5).n,1)];
        end
        MP.BC_psv = [FLD.BC; Nd.BC; false(MP.dof(3).n,1); false(MP.dof(4).n,1);...
        false(MP.dof(5).n,1)];
    end

    MP.tdof = length(MP.BC);
    MP.rdof = sum(MP.BC);

    % for reducing matrix and rhs from reduced active to reduced passive
    MP.act2psv = true(MP.rdof,1);
    MP.act2psv([MP.dof(3).dof_r,MP.dof(4).dof_r,MP.dof(5).dof_r]) = false;
    MP.rdof_psv = sum(MP.act2psv);
end % of function initiate()
