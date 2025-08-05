function MP = set_MP(options)

    global coef

    idir = add_directory_divider(options.input);
    odir = add_directory_divider(options.output);
    MP.idir = idir;
    MP.odir = odir;
    
    if ~exist(odir,'dir') && ~isempty(odir)
        mkdir(odir);
    end

    MP.Estim = 0; % default: 0
    MP.stim_type = options.stim_type;
    % basic parameters taken from the input %
    MP.loc = options.center_location;             % [mm]    
    if options.model_length == 0
        MP.vMC = 1; %
        MP.loc = MP.loc - 1; % 1mm offset to get PoI under slit % [mm]    
        MP.length_BM = 4000; % [um]
    else
        MP.vMC = 0; %
        MP.length_BM = options.model_length*1e3;   % [um]
    end
    MP.extension = 0;       % [um]
    MP.dZ = 10;
    MP.xx = MP.loc-0.5*MP.length_BM*1e-3:0.01:MP.loc+0.5*MP.length_BM*1e-3;
    MP.zz = (-0.5*(MP.length_BM + MP.extension):MP.dZ:0.5*(MP.length_BM + MP.extension))';
    MP.nz = length(MP.zz);
    MP.Nf = options.num_of_freqs;
    if MP.Nf == 0
        MP.freq = options.freq;
        MP.Nf = length(MP.freq);
    else
        MP.freq = transpose(logspace(log10(options.freq(1)),log10(options.freq(2)),MP.Nf));
    end

    MP.BMC = {'AF'};
    MP.Pstim = options.SPL + 25;            % dB SPL at the stapes re. 20 uPa
    if MP.vMC
        MP.Istim = options.stim_current;
    end
    MP.opt_cfld = options.fluid_Corti; % 1: consider corti fluid
    MP.p_alpha = coef.perm_multiplier*options.permeation_coef(1); % permeability constant k/(mu*L) [ms*um^2/ng]
    MP.p_beta = options.permeation_coef(2);
    MP.sensitive_cochlea = options.sensitive_cochlea;
    MP.nonlin = options.nonlin; % 1:nonlinear, 0:linear
    % hard-coded parameters %
    MP.H = 600;                % [um], This may be not used. Instead botH and topH are used.
    MP.CF_H = 50; % Cortfi fluid height [um]
    MP.rho = 1e-3; % fluid density [mg/mm^3]
    MP.nu = coef.nu_multiplier*700; % kinematic viscosity [um^2/ms]
    MP.m = 3; % power of attenuation function to resolve apical wave reflection
    MP.b = 1;

    % model options %
    MP.NQ = 4; % number of gauss points for fluid quadrilateral element
    MP.Visc = 0;        % Indicates viscous (1) or inviscid (0) model
    MP.key_FSI_2Chambers = 3;
    MP.key_FSI = 1;
    MP.opt_wtr_shnt = 0; % 1: consider water shunt
    MP.shnt_r_fact = 1; % water shunt relaxation factor
    MP.ind_MET = 2; % num of MET channel states, either 2 or 10
    MP.opt_in = 1;              % 1: pressure,  2: acceleration
    if options.stim_protocol == 1
        MP.opt_stim = 'multi'; % pure, multi
    elseif options.stim_protocol == 0
        MP.opt_stim = 'pure';
    end
    
    MP.phase = 0.5*pi;              % phase between M-stim and E-stim in degrees     
    
    % sp;ver settings
    MP.gamma = 0.5000;
    MP.beta = 0.2500;
    MP.am = 0.5; 
    MP.af = 0.5;
    MP.cycles = options.cycles;
    MP.T = double(MP.cycles)/min(MP.freq);
    MP.period = 1/min(MP.freq);
    MP.steps = 40;
    MP.save_step = 1; % steps    
    MP.dt = MP.period/MP.steps; % time step
    MP.save_dt = MP.save_step*MP.dt;

    % parameters
    if ~isempty(options.parameter)
        f_names = fieldnames(options.parameter);
        for ii = 1:numel(f_names)
            MP.(f_names{ii}) = options.parameter.(f_names{ii});
        end
    end
end

function pdir = add_directory_divider(pdir)
% % function pdir = add_directory_divider(pdir)
% % If there is no directory divider at the end of the string pdir
% % add it.

    if ~isempty(pdir),
        for ii=1:length(pdir),
            if strcmp(pdir(ii),'\'),
                pdir(ii) = '/';
            end        
        end

        if ~(strcmp(pdir(length(pdir)),'/') || strcmp(pdir(length(pdir)),'\')),
            pdir = [pdir '/'];
        end
    end
    
end