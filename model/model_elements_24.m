function El = model_elements_24(Nd,MP)
% # function El = model_elements(Nd, MP)

    Prop = MP.Prop_fit;
    tilt_OHC = feval(Prop.('tilt_OHC'),MP.xx(:));
    nz = size(Nd.Z,1);
        
    El.Nd1=[]; El.Nd2=[]; El.Nd3=[]; El.Nd4=[]; El.name=[];

    % BM radial
    idx = find(strncmp(Nd.name(1,:),'A',1)); % node names along radial dir.
    [~,s_idx] = sort(Nd.X(1,idx),'ascend');
    idx = idx(s_idx);
    nd1 = (1:nz)' + (idx(1:end-1)-1)*nz;
    nd2 = (1:nz)' + (idx(2:end)-1)*nz;
    El.Nd1 = [El.Nd1; nd1(:)];
    El.Nd2 = [El.Nd2; nd2(:)];
    nAP = find(strcmp(Nd.name,'AP'));
    xE = 0.5*(Nd.X(nd1) + Nd.X(nd2));
    name = repmat({'BMp'},nz*(numel(idx)-1),1);
    name(xE<=Nd.X(nAP)) = {'BMa'};
    El.name = [El.name;name];

    % BM longi.
    nd1 = (1:nz-1)' + (idx-1)*nz;
    nd2 = (2:nz)' + (idx-1)*nz;
    El.Nd1 = [El.Nd1; nd1(:)];
    El.Nd2 = [El.Nd2; nd2(:)];
    nAP = find(strcmp(Nd.name,'AP'));
    xE = 0.5*(Nd.X(nd1) + Nd.X(nd2));
    ap = 0.5*(Nd.X(nAP(1:end-1)) + Nd.X(nAP(2:end)));
    name = repmat({'BMzp'},(nz-1)*numel(idx),1);
    name(xE<=ap) = {'BMza'};
    El.name = [El.name;name];

    % TM radial
    idx = find(strncmp(Nd.name(1,:),'E',1)); % node names along radial dir.
    [~,s_idx] = sort(Nd.X(1,idx),'ascend');
    idx = idx(s_idx);
    nd1 = (1:nz)' + (idx(1:end-1)-1)*nz;
    nd2 = (1:nz)' + (idx(2:end)-1)*nz;
    El.Nd1 = [El.Nd1; nd1(:)];
    El.Nd2 = [El.Nd2; nd2(:)];
    nEA = find(strcmp(Nd.name,'EA'));
    xE = 0.5*(Nd.X(nd1) + Nd.X(nd2));
    name = repmat({'TMx2'},nz*(numel(idx)-1),1);
    name(xE<=Nd.X(nEA)) = {'TMx1'};
    El.name = [El.name;name];

    % TM longi.
    nd1 = (1:nz-1)' + (idx-1)*nz;
    nd2 = (2:nz)' + (idx-1)*nz;
    El.Nd1 = [El.Nd1; nd1(:)];
    El.Nd2 = [El.Nd2; nd2(:)];
    nEA = find(strcmp(Nd.name,'EA'));
    xE = 0.5*(Nd.X(nd1) + Nd.X(nd2));
    ae = 0.5*(Nd.X(nEA(1:end-1)) + Nd.X(nEA(2:end)));
    name = repmat({'TMz1'},(nz-1)*numel(idx),1);
    name(xE<=ae) = {'TMz0'};
    El.name = [El.name;name];

    % OPC radial
    idx = find(strcmp(Nd.name(1,:),'CC') | strncmp(Nd.name(1,:),'CO',2) | ...
               strcmp(Nd.name(1,:),'AP')); % node names along radial dir.
    [~,s_idx] = sort(Nd.Y(1,idx),'ascend');
    idx = idx(s_idx);
    nd1 = (1:nz)' + (idx(1:end-1)-1)*nz;
    nd2 = (1:nz)' + (idx(2:end)-1)*nz;
    El.Nd1 = [El.Nd1; nd1(:)];
    El.Nd2 = [El.Nd2; nd2(:)];
    name = repmat({'OPC'},nz*(numel(idx)-1),1);
    El.name = [El.name;name];

    % IPC radial
    idx = find(strcmp(Nd.name(1,:),'CC') | ...
               strncmp(Nd.name(1,:),'CI',2)); % node names along radial dir.
    [~,s_idx] = sort(Nd.Y(1,idx),'ascend');
    idx = idx(s_idx);
    nd1 = (1:nz)' + (idx(1:end-1)-1)*nz;
    nd2 = (1:nz)' + (idx(2:end)-1)*nz;
    El.Nd1 = [El.Nd1; nd1(:)];
    El.Nd2 = [El.Nd2; nd2(:)];
    name = repmat({'IPC'},nz*(numel(idx)-1),1);
    El.name = [El.name;name];

    % TOC tip longi
    idx = find(strcmp(Nd.name(1,:),'CC'));
    nd1 = (1:nz-1)' + (idx-1)*nz;
    nd2 = (2:nz)' + (idx-1)*nz;
    El.Nd1 = [El.Nd1; nd1(:)];
    El.Nd2 = [El.Nd2; nd2(:)];
    name = repmat({'TCz'},(nz-1)*numel(idx),1);
    El.name = [El.name;name];    

    % RL radial
    idx = find(strcmp(Nd.name(1,:),'CC') | ...
               strncmp(Nd.name(1,:),'B',1)); % node names along radial dir.
    [~,s_idx] = sort(Nd.X(1,idx),'ascend');
    idx = idx(s_idx);
    nd1 = (1:nz)' + (idx(1:end-1)-1)*nz;
    nd2 = (1:nz)' + (idx(2:end)-1)*nz;
    El.Nd1 = [El.Nd1; nd1(:)];
    El.Nd2 = [El.Nd2; nd2(:)];
    name = repmat({'RLx1'},nz*(numel(idx)-1),1);
    El.name = [El.name;name];

    idx = find(strcmp(Nd.name(1,:),'BB') | ...
               strcmp(Nd.name(1,:),'PP')); % node names along radial dir.
    [~,s_idx] = sort(Nd.X(1,idx),'ascend');
    idx = idx(s_idx);
    nd1 = (1:nz)' + (idx(1:end-1)-1)*nz;
    nd2 = (1:nz)' + (idx(2:end)-1)*nz;
    El.Nd1 = [El.Nd1; nd1(:)];
    El.Nd2 = [El.Nd2; nd2(:)];
    name = repmat({'RLp'},nz*(numel(idx)-1),1);
    El.name = [El.name;name];

    % RL longi
    idx = find(strcmp(Nd.name(1,:),'BB'));
    nd1 = (1:nz-1)' + (idx-1)*nz;
    nd2 = (2:nz)' + (idx-1)*nz;
    El.Nd1 = [El.Nd1; nd1(:)];
    El.Nd2 = [El.Nd2; nd2(:)];
    name = repmat({'RLz'},(nz-1)*numel(idx),1);
    El.name = [El.name;name];

    % RLp longi
    idx = find(strcmp(Nd.name(1,:),'PP'));
    nd1 = (1:nz-1)' + (idx-1)*nz;
    nd2 = (2:nz)' + (idx-1)*nz;
    El.Nd1 = [El.Nd1; nd1(:)];
    El.Nd2 = [El.Nd2; nd2(:)];
    name = repmat({'RLpz'},(nz-1)*numel(idx),1);
    El.name = [El.name;name];

    % DCb radial
    idx = find(strcmp(Nd.name(1,:),'AD') | ...
               strcmp(Nd.name(1,:),'DD')); % node names along radial dir.
    [~,s_idx] = sort(Nd.Y(1,idx),'ascend');
    idx = idx(s_idx);
    nd1 = (1:nz)' + (idx(1:end-1)-1)*nz;
    nd2 = (1:nz)' + (idx(2:end)-1)*nz;
    El.Nd1 = [El.Nd1; nd1(:)];
    El.Nd2 = [El.Nd2; nd2(:)];
    name = repmat({'DCb'},nz*(numel(idx)-1),1);
    El.name = [El.name;name];

    % DC longi
    idx = find(strcmp(Nd.name(1,:),'DD'));
    nd1 = (1:nz-1)' + (idx-1)*nz;
    nd2 = (2:nz)' + (idx-1)*nz;
    El.Nd1 = [El.Nd1; nd1(:)];
    El.Nd2 = [El.Nd2; nd2(:)];
    name = repmat({'DCz'},(nz-1)*numel(idx),1);
    El.name = [El.name;name];

    % OHC radial
    idx = find(strcmp(Nd.name(1,:),'DD') | ...
               strcmp(Nd.name(1,:),'BB')); % node names along radial dir.
    [~,s_idx] = sort(Nd.Y(1,idx),'ascend');
    idx = idx(s_idx);
    nd1 = (1:nz)' + (idx(1:end-1)-1)*nz;
    nd2 = (1:nz)' + (idx(2:end)-1)*nz;
    El.Nd1 = [El.Nd1; nd1(:)];
    El.Nd2 = [El.Nd2; nd2(:)];
    name = repmat({'OHC'},nz*(numel(idx)-1),1);
    El.name = [El.name;name];    

    % TOC brackets
    idx = find(strcmp(Nd.name(1,:),'CI1') | ...
               strcmp(Nd.name(1,:),'B1')); % node names along radial dir.
    [~,s_idx] = sort(Nd.Y(1,idx),'ascend');
    idx = idx(s_idx);
    nd1 = (1:nz)' + (idx(1:end-1)-1)*nz;
    nd2 = (1:nz)' + (idx(2:end)-1)*nz;
    El.Nd1 = [El.Nd1; nd1(:)];
    El.Nd2 = [El.Nd2; nd2(:)];
    name = repmat({'IPCb'},nz*(numel(idx)-1),1);
    El.name = [El.name;name];

    idx = find(strcmp(Nd.name(1,:),'CO1') | ...
               strcmp(Nd.name(1,:),'B1')); % node names along radial dir.
    [~,s_idx] = sort(Nd.Y(1,idx),'ascend');
    idx = idx(s_idx);
    nd1 = (1:nz)' + (idx(1:end-1)-1)*nz;
    nd2 = (1:nz)' + (idx(2:end)-1)*nz;
    El.Nd1 = [El.Nd1; nd1(:)];
    El.Nd2 = [El.Nd2; nd2(:)];
    name = repmat({'OPCb'},nz*(numel(idx)-1),1);
    El.name = [El.name;name];

    idx = find(strcmp(Nd.name(1,:),'CO3') | ...
               strcmp(Nd.name(1,:),'A1')); % node names along radial dir.
    [~,s_idx] = sort(Nd.Y(1,idx),'ascend');
    idx = idx(s_idx);
    nd1 = (1:nz)' + (idx(1:end-1)-1)*nz;
    nd2 = (1:nz)' + (idx(2:end)-1)*nz;
    El.Nd1 = [El.Nd1; nd1(:)];
    El.Nd2 = [El.Nd2; nd2(:)];
    name = repmat({'OPCb'},nz*(numel(idx)-1),1);
    El.name = [El.name;name];

    % OTHERS
    ff_r = find(strcmp(Nd.name(1,:),'FF'));
    ex_r = find(strcmp(Nd.name(1,:),'EX'));
    dd_r = find(strcmp(Nd.name(1,:),'DD'));
    bb_r = find(strcmp(Nd.name(1,:),'BB'));
    pp_r = find(strcmp(Nd.name(1,:),'PP'));

    for iz = 1:nz

        ex = (ex_r-1)*nz + iz;
        ff = (ff_r-1)*nz + iz;
        bb = (bb_r-1)*nz + iz;
        pp = (pp_r-1)*nz + iz;
        dd = (dd_r-1)*nz + iz;

        bbp = pp + tilt_OHC(iz); % bbp: tilt of DCp toward the apical direction by distance of tilt_OHC
        bbq = pp - tilt_OHC(iz); % bbp: tilt of DCp toward the basal direction by distance of tilt_OHC

        if tilt_OHC(iz)>=0 && nz-iz >= tilt_OHC(iz)
            El.Nd1 = [El.Nd1; dd];
            El.Nd2 = [El.Nd2; bbp];
            El.name = [El.name; {'DCp'}];         
        elseif tilt_OHC(iz)<=0 && iz >= -tilt_OHC(iz)                
            El.Nd1 = [El.Nd1; dd];
            El.Nd2 = [El.Nd2; bbp];
            El.name = [El.name; {'DCp'}]; 
        end            
                    
        if tilt_OHC(iz)>=0 && nz-iz < tilt_OHC(iz)   % end of triangular repeats, suppressed
            El.Nd1 = [El.Nd1; dd];
            El.Nd2 = [El.Nd2; bbq];
            El.name = [El.name; {'DCp'};]; % end of triangular repeats, suppressed
        elseif tilt_OHC(iz)<=0 && iz < tilt_OHC(iz)
            El.Nd1 = [El.Nd1; dd;];
            El.Nd2 = [El.Nd2; bbq;];
            El.name = [El.name; {'DCp'}];            
        end

        % % OHC bundle and linear spring
        El.Nd1 = [El.Nd1; bb;];
        El.Nd2 = [El.Nd2; ex;];
        El.name = [El.name; {'OHB'};];

        % % Pseudo rod at OHC bundle location
        El.Nd1 = [El.Nd1; bb;];
        El.Nd2 = [El.Nd2; ff;];
        El.name = [El.name; {'xHB'};];
        
        % % Pseudo link to conveniantly calculate translational displ. 
        % % or damping of HB            
        El.Nd1 = [El.Nd1; ff;];
        El.Nd2 = [El.Nd2; ex;];
        El.name = [El.name; {'ANK'};]; 
    end

    El.N = length(El.Nd1);  

end