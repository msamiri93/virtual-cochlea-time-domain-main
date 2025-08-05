function CFLD = model_CFLD(MP,Nd)

    CFLD.L = MP.length_BM;
    CFLD.H = MP.CF_H;

    CFLD.rho = MP.rho;                     % fluid density [mg/mm^3]

    zz = (-0.5*MP.length_BM:MP.dZ:0.5*MP.length_BM)';
    idx = knnsearch(MP.zz(:),zz(:));

%     CFLD.nBot = find(strcmp(Nd.name,'AM')); % BM centerline nodes number
%     CFLD.nBot = find(strcmp(Nd.name,MP.BMC)); % BM centerline nodes number
    CFLD.nBot = zeros(MP.nz,numel(MP.BMC));
    for i = 1:numel(MP.BMC)
        CFLD.nBot(:,i) = find(strcmp(Nd.name,MP.BMC{i})); % BM centerline nodes number
    end    
    CFLD.nTop = find(strcmp(Nd.name,'BB')); %
%     CFLD.nTop = [find(strcmp(Nd.name,'BB')),find(strcmp(Nd.name,'CC'))]; %
%     CFLD.nTop = find(strcmp(Nd.name,'CC')); % 
%     CFLD.nTop = find(strcmp(Nd.name,'E1'));
    CFLD.nBot = CFLD.nBot(idx,:); 
    CFLD.nTop = CFLD.nTop(idx,:);
    CFLD.Nd.width = Nd.X(CFLD.nBot);
%     CFLD.dcf_32 = 0.6;                   % Dimension conversion factor from 3-D to 2-D (acceleration felt by fluid boudnary)
%     CFLD.dcf_23 = 0.6;                   % Dimension conversion factor from 2-D to 3-D (pressure felt by CP strutural surfaces)
    CFLD.idx = idx;
end