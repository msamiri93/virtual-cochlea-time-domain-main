function [O,Nd,El,AMat,MP,FLD,OHC,C,c] = virtualCochleaDynamic(options)
% [R,Nd,El,AMat,MP,FLD,OHC,U,uu] = virtualCochleaDynamic(options)
%
% Arguments:
%    input                 Input directory
%    output                Output directory
%    center_location       Center location in [mm]
%    model_length          L of modeled piece in [mm] - model_length = 0: micro-chamber model
%    freq                  Stim freq range [f1, f2], single frequency, or multiple frequencies stimulation in [kHz]
%    num_of_freqs          Number of frequencies [0: only input frequencies]
%    SPL                   Stim pressure level in [dB SPL]
%    permeation_coef       Permeation coef [Alpha, Beta]
%    stim_protocol         Stimulation protocol [0:Pure tone, 1:Multi tone]
%    fluid_Corti           Corti fluid [1:On, 0:Off]
%    sensitive_cochlea     Sensitive cochlea [1:On, 0:Off]
%    nonlin                Nonlinearity [1:On, 0:Off]
%    stim_current          Stim current level in [mA]     
%    stim_type
%
% Example:
% [R,Nd,El,MP,FLD,OHC,U,uu] = virtualCochleaDynamic(freq = 4, fluid_Corti = 1, sensitive_cochlea = 1, nonlin = 1, SPL = 40, cycles = 20, plot = true, save = true);

    arguments
        options.input (1,:) char = ' ';                     % Input directory
        options.output (1,:) char = './houtput';            % Output directory
        options.center_location (1,1) double       = 6;     % Center location in [mm]
        options.model_length (1,1) double          = 12;    % L of modeled piece in [mm]
        options.freq (1,:) double       = 1;                % Stim freq range [f1, f2] in [kHz]
        options.num_of_freqs (1,1) int16 = 0;               % Num of frequencies [0: only input]
        options.cycles (1,1) int16 = 15;                    % Num of cycles
        options.SPL (1,1) double = 40;                      % Stim pressure level in [dB SPL]
        options.permeation_coef (1,2) double = [1e3,1];     % Permeation coef [Alpha, Beta]
        options.stim_protocol (1,1) int16 = 0;              % Stimulation protocol [0:Pure tone, 1:Multi tone]
        options.fluid_Corti (1,1) int16 = 1;                % Corti fluid [1:On, 0:Off]
        options.sensitive_cochlea (1,1) int16 = 1;          % Sensitive cochlea
        options.nonlin (1,1) int16 = 0;                     % Nonlinearity
        options.stim_current (1,1) double = 0.1;            % Stim current level in [mA]     
        options.stim_type {mustBeMember(options.stim_type,{'impulse','tone','random'})} = 'random';
        options.log = false;
        options.plot = false;
        options.save = false;

        options.parameter = [];
    end

    set_Global_coefficients; % this is setting coefficients for parameter study
    model_coefs;
    % add_paths;

    To = onset_msg(options.log);

    global coef  

    MP = set_MP(options);
    MP.coef = coef; % save coef information before save

    [Nd,El,OHC,IHC,MP,FLD] = initiate(MP);

    % frequency independent part. only contains inviscid scala fluid formulation
    [A0,R0,~,R,AMat,El,FLD,OHC,~,Ge,Ie,V0] = assemble_A0(Nd,El,OHC,IHC,FLD,MP);

    Vnd = reshape([OHC.Vnd],5,[]); Vm0 = Vnd(2,:)-Vnd(3,:);
    HB = [OHC.HB]; po0 = [HB.po0];

    dA = decomposition(A0);
    clear A0;

    U = initiate_U(MP,OHC);
    an = U(MP.dof(2).udof,1);
    vn = U(MP.dof(2).udof,1);
    xn = U(MP.dof(2).udof,1);
    uu = zeros(MP.dof(2).n,MP.steps + 1);

    xx_idx = 1:4:MP.nz;
    Z0 = NaN(numel(xx_idx),MP.steps*MP.cycles/MP.save_step + 1);
    O = struct('yBM_t',Z0,'yRL_t',Z0,'yDC_t',Z0,'xDC_t',Z0,'xHB_t',Z0,...
               'po_t',Z0,'Vm_t',Z0,'time',NaN(1,MP.steps*MP.cycles/MP.save_step + 1));
    O.xx = MP.xx(xx_idx);
    clear Z0;

    b0 = initiate_RHS(MP,Nd,OHC,AMat,Ge,V0,Ie); % initiate the vector of variables and RHS vector

    time = 0;
    icount = 0;
    tcount = 0;
    peak1 = eps; peak2 = eps; peak3 = eps;

    nAF = find(strcmp(Nd.name,'AF'));
    nBB = find(strcmp(Nd.name,'BB'));
    nDD = find(strcmp(Nd.name,'DD'));
    nANK = find(strcmp(El.name,'ANK'));
    adof = MP.dof(4).adof;
    odof = MP.dof(5).odof;
    edof = reshape(MP.dof(3).edof,5,MP.nz);
    nd1 = El.Nd1(nANK); nd2 = El.Nd2(nANK);

    f_names = fieldnames(O);
    for ii = 1:numel(f_names);
        switch f_names{ii}
            case 'po_t'
                O.(f_names{ii})(:,1) = U(odof(xx_idx),1);
            case 'Vm_t'
                O.(f_names{ii})(:,1) = U(edof(2,xx_idx),1) - U(edof(3,xx_idx),1);
            otherwise
                O.(f_names{ii})(:,1) = 0;
        end   
    end

    if options.plot
        % probing response
        figure; set(gcf,'Position',1e3*[0.6,0.5,0.84,0.4]);
        t = tiledlayout(gcf,2,3);
        ax = gobjects(2,3);
        ax(1,1) = nexttile(t); ax(1,2) = nexttile(t); ax(1,3) = nexttile(t);
        ax(2,1) = nexttile(t); ax(2,2) = nexttile(t); ax(2,3) = nexttile(t);
        % probing input/output
        figure; set(gcf,'Position',1e3*[0.62,0.52,0.84,0.24]);
        t = tiledlayout(gcf,1,3);
        io_ax = gobjects(1,3);
        io_ax(1,1) = nexttile(t); io_ax(1,2) = nexttile(t); io_ax(1,3) = nexttile(t);
        h_input = animatedline(io_ax(1,1));
        h_output = animatedline(io_ax(1,2));
        h_err = animatedline(io_ax(1,3)); yscale(io_ax(1,3),'log');
    end

    if options.save
        ct = clock;
        save_data(MP,'MP',Nd,'Nd',El,'El',OHC,'OHC',FLD,'FLD', ...
            MP = MP, name = 'model', output = options.output, SPL = options.SPL, ct = ct);
        % save(file_name,'MP','Nd','El','OHC','FLD','AMat');
    end

    dt_s = 0.1*MP.dt;
    % Implicit channel kinetics
    dxadt = zeros(MP.dof(4).n,1);
    dpodt = zeros(MP.dof(4).n,1);
    %
    interruptFlag = false; % Flag to check for interruption
    %
    convergence_flag = false;
    rel_tol = 1e-2;
    %
    while time < MP.T
        try
            icount = icount + 1;
    
            a = U(MP.dof(2).udof,icount);        
            xn = xn + MP.dt*vn + (MP.dt^2)*((0.5-MP.beta)*an + MP.beta*a);
            vn = vn + MP.dt*((1-MP.gamma)*a + MP.gamma*an);
            an = a;
            
            xx = reshape(xn,6,Nd.N); yBM = xx(2,nAF);
            xHB = xx(1,nd2) - xx(1,nd1);
            xa = U(adof,icount);
            po = U(odof,icount);
            dpo = po - po0.';
            Vnd = U(edof,icount);
            Vm = U(edof(2,:),icount)-U(edof(3,:),icount);
    
            if options.plot && (mod(time,4*MP.dt) < MP.dt)
                if 1.1*max(abs(yBM)) > peak1, peak1 = 1.1*max(abs(yBM)); end
                cla(ax(1,1)); plot(ax(1,1),MP.xx,yBM); ylim(ax(1,1),[-peak1,peak1]);
                hold(ax(1,1),"on"); plot(ax(1,1),MP.xx,xHB); plot(ax(1,1),MP.xx,xa);
        
                if 1.1*max(abs(po)) > peak2, peak2 = 1.1*max(abs(po)); end
                cla(ax(1,2)); plot(ax(1,2),MP.xx,po); ylim(ax(1,2),[0,1]);
        
                if 1.1*max(abs(Vm-Vm0.')) > peak3, peak3 = 1.1*max(abs(Vm-Vm0.')); end
                cla(ax(1,3)); plot(ax(1,3),MP.xx,Vm-Vm0.'); ylim(ax(1,3),[-peak3,peak3]);
                cla(ax(2,1)); plot(ax(2,1),O.time,O.yBM_t(knnsearch(O.xx',2),:),O.time,O.yRL_t(knnsearch(O.xx',2),:));
                cla(ax(2,2)); plot(ax(2,2),O.time,O.yBM_t(knnsearch(O.xx',6),:),O.time,O.yRL_t(knnsearch(O.xx',6),:));
                cla(ax(2,3)); plot(ax(2,3),O.time,O.yBM_t(knnsearch(O.xx',10),:),O.time,O.yRL_t(knnsearch(O.xx',10),:));                   
            end
    
            [b, pre_in] = evaluate_RHS(time,b0,MP,Nd,El,FLD,OHC,AMat,R0,R,U(:,icount),an,vn,xn,dpo,Vnd);
    
            u = dA\b(MP.BC,1);
    
            U(:,icount+1) = update_U(time,MP,Nd,FLD,u,U(:,icount+1), pre_in);
    
            if MP.nonlin
                % [xa, po, dt_s] = explicit_channel_kinetics_r02(time, dt_s, MP, HB, po, xHB, xa);
                [xa, po, dxadt, dpodt, dt_s] = implicit_channel_kinetics(time, dt_s, MP, HB, po, xHB, xa, dxadt, dpodt);
                U(adof,icount+1) = xa;
                U(odof,icount+1) = po;
            end
            
            a = U(MP.dof(2).udof,icount+1);
            x = xn + MP.dt*vn + (MP.dt^2)*((0.5-MP.beta)*an + MP.beta*a);
    
            uu(:,icount+1) = x;
    
            time = time + MP.dt;
    
            xx = reshape(uu(:,icount+1),6,Nd.N);
    
            tcount = tcount + 1;
            O.yBM_t(:,tcount + 1) = xx(2,nAF(xx_idx));
            O.yRL_t(:,tcount + 1) = xx(2,nBB(xx_idx));
            O.xDC_t(:,tcount + 1) = xx(1,nDD(xx_idx));
            O.yDC_t(:,tcount + 1) = xx(2,nDD(xx_idx));
            O.xHB_t(:,tcount + 1) = xx(1,nd2(xx_idx)) - xx(1,nd1(xx_idx));
            O.po_t(:,tcount + 1) = U(odof(xx_idx),icount+1);
            O.Vm_t(:,tcount + 1) = U(edof(2,xx_idx),icount+1) - U(edof(3,xx_idx),icount+1);
            O.time(tcount + 1) = time;    

            if options.plot
                if ~exist('sum_of_sq_h_in','var')
                    sum_of_sq_h_in = 0;
                    sum_of_sq_h_out = 0;
                    ref = NaN;
                end
                sum_of_sq_h_in = sum_of_sq_h_in + sum(pre_in.^2);
                sum_of_sq_h_out = sum_of_sq_h_out + sum(xx(2,nAF(xx_idx)).^2);                
                if eq(icount,MP.steps) 
                    addpoints(h_input,time,sqrt(sum_of_sq_h_in/tcount));
                    addpoints(h_output,time,sqrt(sum_of_sq_h_out/tcount));
                    err = abs(sqrt(sum_of_sq_h_out/tcount) - ref)/ref;
                    addpoints(h_err,time,err);
                    ref = sqrt(sum_of_sq_h_out/tcount);
                    if err < rel_tol, convergence_flag = true; end
                end
            end

            if  (icount+1) > MP.steps
                icount = 0;
                C = U;
                c = uu;
                U(:,1) = U(:,end);
                uu(:,1) = x;
                if convergence_flag
                    break;
                end
            end

            if options.plot && (mod(time,4*MP.dt) < MP.dt)
                drawnow
            end

            if interruptFlag
                fprintf('Saving and exiting...\n');
                break; % Exit loop
            end
    
        catch ME
            % --- Handle Errors or Interrupts Gracefully ---
            if strcmp(ME.identifier, 'MATLAB:operandsInterrupted')
                fprintf('Interrupt detected. Saving and exiting...\n');
                interruptFlag = true;
                continue; % Skip to next iteration or break
            else
                rethrow(ME); % Genuine error (not an interrupt)
            end
        end
    end

    if options.save
        cycle = round(time/MP.period);
        cycle_name = ['cycle_',num2str(cycle,'%d')];
        if ~exist('C','var'), C = []; c = []; end
        save_data(MP,'MP',O,'O',C,'C',c,'c', MP = MP, name = cycle_name, ...
            output = options.output, SPL = options.SPL, ct = ct);
    end

    fin_msg(To,options.log);

end

function To = onset_msg(log)
    To=clock;
    if ~log
        return
    end
    fprintf(1,'\n\n   Performing organ of Corti FE Harmonic Analysis...\n\n');
    fprintf(1,'  Program starts at %d:%d\n\n',To(4:5));
end % of function onset_msg()

function fin_msg(To,log)
    if ~log
        return
    end
    Tf = clock;
    eT=etime(Tf,To);
    secs=rem(round(eT),60);
    mins=rem(floor(eT/60),60);
    hrs=floor(eT/3600);
    
    fprintf(1,'  Program ends at %3d:%3d:%3d sec \n',round(Tf(4:6)));
    fprintf(1,'  It took %d hr %d min %d sec total.',hrs,mins,secs);
    fprintf(1,'\n=======================================================\n');
end % of function fin_msg()

function save_data(data,tag,options)

    arguments (Repeating)
        data
        tag
    end

    arguments
        options.MP
        options.SPL
        options.output
        options.name
        options.ct
    end

    MP = options.MP;
        
    ct = options.ct;
    date_str = [num2str(ct(2),'%02.2d'),num2str(ct(3),'%02.2d'),num2str(rem(ct(1),2000),'%02.2d')];
    time_str = [num2str(ct(4),'%02.2d'),num2str(ct(5),'%02.2d')];

    if MP.sensitive_cochlea, status = 'act'; else, status = 'psv'; end

    if numel(MP.freq) > 1
        fdir = [options.output,'/r',date_str,'_',time_str,'_',status];
    else
        freq_str = [num2str(MP.freq,'%02.2d'),'kHz'];
        SPL_str = [num2str(options.SPL,'%02.2d'),'dB'];
        fdir = [options.output,'/r',date_str,'_',time_str,'_',freq_str,'_',SPL_str,'_',status];

    end
    if ~exist(fdir,'dir')
        mkdir(fdir);
    end
    file_name = [fdir,'/',options.name,'.mat'];

    S = cell2struct(data',tag');
    save(file_name,'-struct','S')
end