function [MP,Nd] = model_nodes_24(MP, opt_lgradient)
% # function [MP,Nd] = model_nodes_24(MP)

    Prop_fit = MP.Prop_fit;
    
    p = define_key_coordinates(MP,Prop_fit, opt_lgradient);

    Nd = assign_nodal_coordinates(MP, Prop_fit, p);    

    Nd = initialize_nodes(Nd);

    MP.Prop_fit = Prop_fit;

end

function Nd = assign_nodal_coordinates(MP, Prop, p)

    p_names = fieldnames(p);

    if isfield(MP,'nz')
        nz = MP.nz;
    else
        nz = length(MP.xx);
    end

    xx = MP.xx(:);

    % extra nodes to capture higher order vibrations
    nBMP = 1;
    nOPC = 3;
    nTMB = 1;
    nRL = 1; % changing this value requires changing the code
    nIPC = 1; % changing this value requires changing the code
    nBMA = 1; % changing this value requires changing the code
    nAF = 1; % changing this value requires changing the code
    point_num = numel(p_names);
    node_num = point_num + nAF + nBMA + nBMP + nOPC + nIPC + nTMB + nRL;
  
    % Nd structure dimensions
    dim = {'X','Y','Z'};

    for i = 1:numel(dim)
        Nd.(dim{i}) = zeros(nz,node_num);
    end
    
    Nd.name = cell(nz,node_num);

    for iz = 1:nz
        for pi = 1:point_num
            for di = 1:2
                name = p_names{pi};
                Nd.(dim{di})(iz,pi) = p.(name)(iz,di);
            end
            Nd.name(iz,pi) = {name};
        end
    end

    % z_tilt = -0.25*feval(Prop.('tilt_OHC'),xx)*MP.dZ;
    z_tilt = -0.125*feval(Prop.('tilt_OHC'),xx)*MP.dZ;      % O--+------X

    if MP.dim == 2
        z_tilt = 0;
    end

    idx = strncmp(Nd.name,'E',1) | strncmp(Nd.name,'BB',2) | ...
          strncmp(Nd.name,'P',1) | strncmp(Nd.name,'F',1);
    for iz = 1:nz
        Nd.Z(iz,idx(1,:)) = z_tilt(iz);
    end
    % discretization

    for iz = 1:nz

        cnt = point_num;
        BMi = 1;
        % BMA
        for di = 1:2
            % dd = linspace(p.A0(iz,di),p.AP(iz,di),nBMA+2); dd([1,end]) = [];
            dd = 0.25*p.A0(iz,di) + 0.75*p.AP(iz,di);
            Nd.(dim{di})(iz,cnt+1:cnt+nBMA) = dd;
        end
        idx = strcmp(Nd.name(iz,:),'A0') | strcmp(Nd.name(iz,:),'AP');
        dd = Nd.Z(iz,idx);
        dd = linspace(dd(1),dd(2),nBMA+2); dd([1,end]) = [];
        Nd.Z(iz,cnt+1:cnt+nBMA) = dd;
        names = cell(nBMA,1);
        for i = 1:nBMA
            names{i} = ['A',num2str(BMi)];
            BMi = BMi + 1;
        end
        Nd.name(iz,cnt+1:cnt+nBMA) = names; 
        cnt = cnt + nBMA;
        % BM AP-AD mid node
        for di = 1:2
            dd = linspace(p.AP(iz,di),p.AD(iz,di),3); dd([1,end]) = [];
            Nd.(dim{di})(iz,cnt+1:cnt+1) = dd;
        end
        idx = strcmp(Nd.name(iz,:),'AP') | strcmp(Nd.name(iz,:),'AD');
        dd = Nd.Z(iz,idx);
        dd = linspace(dd(1),dd(2),3); dd([1,end]) = [];
        Nd.Z(iz,cnt+1:cnt+1) = dd;
        names = cell(1,1);
        names{i} = 'AF';
        Nd.name(iz,cnt+1:cnt+1) = names;
        cnt = cnt + 1;
        % BMp
        for di = 1:2
            dd = linspace(p.AD(iz,di),p.AX(iz,di),nBMP+2); dd([1,end]) = [];
            Nd.(dim{di})(iz,cnt+1:cnt+nBMP) = dd;
        end
        idx = strcmp(Nd.name(iz,:),'AD') | strcmp(Nd.name(iz,:),'AX');
        dd = Nd.Z(iz,idx);
        dd = linspace(dd(1),dd(2),nBMP+2); dd([1,end]) = [];
        Nd.Z(iz,cnt+1:cnt+nBMP) = dd;
        names = cell(nBMP,1);
        for i = 1:nBMP
            names{i} = ['A',num2str(BMi)];
            BMi = BMi + 1;
        end
        Nd.name(iz,cnt+1:cnt+nBMP) = names;
        cnt = cnt + nBMP;
        % OPC
        for di = 1:2
            dd = linspace(p.CC(iz,di),p.AP(iz,di),nOPC+2); dd([1,end]) = [];
            Nd.(dim{di})(iz,cnt+1:cnt+nOPC) = dd; 
        end
        idx = strcmp(Nd.name(iz,:),'CC') | strcmp(Nd.name(iz,:),'AP');
        dd = Nd.Z(iz,idx);
        dd = linspace(dd(1),dd(2),nOPC+2); dd([1,end]) = [];
        Nd.Z(iz,cnt+1:cnt+nOPC) = dd;
        names = cell(nOPC,1);
        for i = 1:nOPC
            names{i} = ['CO',num2str(i)];
        end
        Nd.name(iz,cnt+1:cnt+nOPC) = names;
        cnt = cnt + nOPC;
        % TMB
        for di = 1:2
            dd = linspace(p.EA(iz,di),p.EX(iz,di),nTMB+2); dd([1,end]) = [];
            Nd.(dim{di})(iz,cnt+1:cnt+nTMB) = dd;
        end
        idx = strcmp(Nd.name(iz,:),'EA') | strcmp(Nd.name(iz,:),'EX');
        dd = Nd.Z(iz,idx);
        dd = linspace(dd(1),dd(2),nTMB+2); dd([1,end]) = [];
        Nd.Z(iz,cnt+1:cnt+nTMB) = dd;
        names = cell(nTMB,1);
        for i = 1:nTMB
            names{i} = ['EB',num2str(i)];
        end
        Nd.name(iz,cnt+1:cnt+nTMB) = names;
        cnt = cnt + nTMB;
        % RL
        for di = 1:2
            % dd = linspace(p.CC(iz,di),p.BB(iz,di),nRL+2); dd([1,end]) = [];
            dd = 0.35*p.CC(iz,di) + 0.65*p.BB(iz,di);
            Nd.(dim{di})(iz,cnt+1:cnt+nRL) = dd;
        end
        idx = strcmp(Nd.name(iz,:),'CC') | strcmp(Nd.name(iz,:),'BB');
        dd = Nd.Z(iz,idx);
        dd = linspace(dd(1),dd(2),nRL+2); dd([1,end]) = [];
        Nd.Z(iz,cnt+1:cnt+nRL) = dd;
        names = cell(nRL,1);
        for i = 1:nRL
            names{i} = ['B',num2str(i)];
        end
        Nd.name(iz,cnt+1:cnt+nRL) = names;
        cnt = cnt + nRL;
        % IPC      
        for di = 1:2
            dd = linspace(p.CC(iz,di),p.CI(iz,di),nIPC+2); dd([1,end]) = [];
            Nd.(dim{di})(iz,cnt+1:cnt+nIPC) = dd;
        end
        idx = strcmp(Nd.name(iz,:),'CC') | strcmp(Nd.name(iz,:),'CI');
        dd = Nd.Z(iz,idx);
        dd = linspace(dd(1),dd(2),nIPC+2); dd([1,end]) = [];
        Nd.Z(iz,cnt+1:cnt+nIPC) = dd;
        names = cell(nIPC,1);
        for i = 1:nIPC
            names{i} = ['CI',num2str(i)];
        end
        Nd.name(iz,cnt+1:cnt+nIPC) = names;        
    end

    Nd.Z = Nd.Z + xx*1e3;

end

function p = define_key_coordinates(MP,Prop,opt_lgradient)

    if isfield(MP,'nz')
        nz = MP.nz;
    else
        nz = length(MP.xx);
    end

    if opt_lgradient
        xx = MP.xx(:); % in mm
    else
        xx = MP.loc*ones(nz,1);
    end
    
    if sum(feval(Prop.('alpha_TM'),xx) > feval(Prop.('alpha_RL'),xx)) > 0
        error('TM and RL relative angle not acceptable.')
    end

    p.A0 = zeros(nz,2); % root of IPC
    p.AP = zeros(nz,2);
    p.AP(:,1) = p.A0(:,1) + feval(Prop.('width_BMA'),xx);

    p.AD = zeros(nz,2);
    % p.AD(:,1) = p.AP(:,1) + feval(Prop.('r_DC_root'),xx).*feval(Prop.('width_OHC_RL'),xx);
    p.AD(:,1) = p.AP(:,1) + 1.1*feval(Prop.('width_OHC_RL'),xx); % 1 to 1.1 RL and DC root radial extent

    % p.AX = zeros(nz,2);
    % p.AX(:,1) = p.AD(:,1) + feval(Prop.('DC_to_BMLateral'),xx);

    p.AX = zeros(nz,2);
    p.AX(:,1) = p.AP(:,1) + feval(Prop.('width_BMP'),xx);    

    p.CC = zeros(nz,2);
    p.CC(:,2) = feval(Prop.('height_TOC'),xx);
    x_offset = p.CC(:,2)./tand(feval(Prop.('alpha_TOC'),xx));
    p.CC(:,1) = feval(Prop.('IPC_root_offset'),xx) + x_offset;
    p.CI(:,1) = p.CC(:,1) - feval(Prop.('r_IPC'),xx).*x_offset;
    p.CI(:,2) = p.CC(:,2) - feval(Prop.('r_IPC'),xx).*feval(Prop.('height_TOC'),xx);

    a_RL = feval(Prop.('alpha_RL'),xx);
    p.BB = p.CC + [cosd(a_RL),sind(a_RL)].*feval(Prop.('width_OHC_RL'),xx);

    % p.PP = p.CC + [cosd(a_RL),sind(a_RL)].*(feval(Prop.('width_OHC_RL'),xx) + feval(Prop.('diam_DCpr'),xx)); % idk why diam_DCpr was used here
    
    p.PP = p.CC + [cosd(a_RL),sind(a_RL)].*(1.2*feval(Prop.('width_OHC_RL'),xx));

    OHC_DC = vecnorm([p.AD(:,2) - p.BB(:,2),p.AD(:,1) - p.BB(:,1)],2,2);
    a = feval(Prop.('DD_theta'),xx) + atan2d(p.AD(:,2) - p.BB(:,2),...
        p.AD(:,1) - p.BB(:,1));

    p.DD = p.BB + 0.4*[cosd(a), sind(a)].*OHC_DC; % uniform 0.4 OHC, 0.6 DC

    a = 90 + a_RL; % hair bundle tilt can be added later, for now it is perpendicular to RL 
    p.EX = p.BB + [cosd(a - 1), sind(a - 1)].*feval(Prop.('height_HB'),xx);
    p.FF = p.BB + [cosd(a + 4), sind(a + 4)].*feval(Prop.('height_HB'),xx);
    a = 180 + feval(Prop.('alpha_TM'),xx);
    p.E0 = p.EX + [cosd(a), sind(a)].*feval(Prop.('width_TM'),xx);
    p.EA = (3/4)*p.E0 + (1/4)*p.EX;

end

function Nd = initialize_nodes(Nd)
        
    Nd.N = length(Nd.X);
    Nd.Nc = size(Nd.X,2);
    Nd.Nr = size(Nd.X,1);
        
    Nd.dx = zeros(size(Nd.X));
    Nd.dy = zeros(size(Nd.X));
    Nd.dz = zeros(size(Nd.X));
    Nd.rx = zeros(size(Nd.X));
    Nd.ry = zeros(size(Nd.X));
    Nd.rz = zeros(size(Nd.X));
    
    NDOF = 6;
    Nd.BC = ones(Nd.Nr* Nd.Nc, NDOF);
    Nd.N = Nd.Nr*Nd.Nc;
end