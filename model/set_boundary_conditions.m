function Nd = set_boundary_conditions(MP, Nd, opt_MP_set)

    if ~exist('opt_MP_set','var') %  a new set of geometry and material props updated on Apr 2021
        opt_MP_set = 1;         % default is using the new set. 0 = the set before 2021
    end

    nd_SLam = find(strcmp(Nd.name,'A0'));       % Spiral laminar at IPC foot (medial edge)
    nd_SLig = find(strcmp(Nd.name,'AX'));       % Spiral ligament (lateral edge)
    nd_TLam = find(strcmp(Nd.name,'E0'));       % Spiral lamina at TM attach't    
    
    nd_IPCToe = find(strcmp(Nd.name,'CI'));     % IPC toe
%     nd_IPCToe = find(strcmp(Nd.name,'A2'));     % IPC toe

    [nr,nc] = size(Nd.X);    
%     nd_bend = (1:nr:nr*nc);
    locz = -0.5*(MP.length_BM); % [um]
    nd_bend = find(Nd.Z<0.4*MP.dZ + locz & Nd.Z>=-0.5*MP.dZ + locz);
    locz = 0.5*(MP.length_BM); % [um]
    nd_aend = find(Nd.Z<0.4*MP.dZ + locz & Nd.Z>=-0.5*MP.dZ + locz);

%     nd_aend = (nr:nr:nr*nc);
    %nd_aend = [(nr-1:nr:(nr-1)*nc), (nr:nr:nr*nc)];
    % if MP.vMC == 1
    %     nd_ends = [nd_bend; nd_aend];
    % else
        nd_ends = [nd_bend;nd_aend];
    % end
    if opt_MP_set == 0    
        Nd = set_BC3(Nd, nd_SLam, 'dAll');
    else
        Nd = set_BC3(Nd, nd_SLam, 'All');        
        Nd = set_BC3(Nd, nd_IPCToe, 'All');
        Nd = set_BC3(Nd, nd_IPCToe, 'rx');
        Nd = set_BC3(Nd, nd_IPCToe, 'ry');        
    end
    Nd = set_BC3(Nd, nd_SLig, 'All');
    
    Nd = set_BC3(Nd, nd_TLam, 'All');

    Nd = set_BC3(Nd, nd_ends, 'dAll'); 
    % Nd = set_BC3(Nd, nd_ends, 'dz');    
%     Nd = set_BC3(Nd, nd_ends, 'dx');
    % Nd = set_BC3(Nd, nd_ends, 'rx');
    % Nd = set_BC3(Nd, nd_ends, 'ry');
    % Nd = set_BC3(Nd, nd_ends, 'rz');  

    % hinge condition
    Nd = set_BC3(Nd,Nd.Master_Slave(:,2),'dAll');
    
    Nd = find_mapF2R2(Nd);
    
    % MSA: this was adjusted in the end to make using Nd.BC easier 
    Nd.BC = reshape(transpose(logical(Nd.BC)),Nd.tdof,1);
    
end



function Nd = set_BC3(Nd, cnode, cdof)

% %      Nd = set_boundary_condition(Nd, cnode, cdof)
% % 
% %      Description: assign boundary conditions
% % 
% %      Nd: node structure, cnode: node id to constrain, 
% %      cdof: DOF to confine (1=dx, 4=rx, ...)
% %     
% %      First written: 2009.03.27
% %      Last modified: 
% %      by Jong-Hoon Nam
% % 

    if ischar(cnode)
        if strcmpi(cnode(1),'A')
            cnode = 1:Nd.N;
        end
    end
    
    if ~isempty(cnode)
    
        switch (upper(cdof(1)))
            case 'A'
                Nd.BC(cnode,:) = 0;
            case 'D'
                if strcmpi(cdof(2),'X')
                    Nd.BC(cnode,1) = 0;
                elseif strcmpi(cdof(2),'Y')
                    Nd.BC(cnode,2) = 0;
                elseif strcmpi(cdof(2),'Z')
                    Nd.BC(cnode,3) = 0;
                elseif strcmpi(cdof(2),'A')
                    Nd.BC(cnode,1:3) = 0;
                else
                    printf(1,'\n Warning: Unknown option for displ constraint \n');
                end
            case 'R'
                if strcmpi(cdof(2),'X')
                    Nd.BC(cnode,4) = 0;
                elseif strcmpi(cdof(2),'Y')
                    Nd.BC(cnode,5) = 0;
                elseif strcmpi(cdof(2),'Z')
                    Nd.BC(cnode,6) = 0;
                elseif strcmpi(cdof(2),'A')
                    Nd.BC(cnode,4:6) = 0;
                else
                    printf(1,'\n Warning: Unknown option for rotational constraint \n');
                end            
            otherwise
                fprintf(1,'\n Warning: check if the boundary condition was properly assigned. \n');
        end
        
    end
end