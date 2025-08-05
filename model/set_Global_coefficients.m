function set_Global_coefficients(z)
% function set_Global_coefficients(z)
% z: structure containing parameters; usually used for reruns
%
    
    clear global;
    global coef

    if exist('z','var')
        coef = z;
    else
        % gating swing
        coef.HB_b_A = 1.3; % [nm]
        coef.HB_b_B = 1.3; % [nm]
        
        %fluid
        coef.dcf_32_A = 0.6;
        coef.dcf_23_A = 0.6;
    
        coef.dcf_32_B = 0.6;
        coef.dcf_23_B = 0.6;
    
        %damping
        d_opt = 0;
        coef.beta_multiplier = 1;
        coef.nu_multiplier = 10; % Corti fluid nu
        if d_opt == 1
            coef.BM_damping_coef = 0.9;
            coef.RL_damping_coef = 0.9;
            coef.TM_damping_coef = 0.9;
            coef.DC_damping_coef = 0.9;
            coef.PC_damping_coef = 0.9;
        end   
        %properties
        coef.rho_multiplier = 1; % fluid density
        coef.perm_multiplier = 0.1;
    end
        coef.hinge_ratio = 0.9;
        coef.OHC_damping = 1;
end