function FLD = model_FLD_cochlea(MP,H,Nd)                     % size of helicotrema opening

    FLD.L = MP.length_BM;
    FLD.H = H;

    FLD.rho = MP.rho;                     % fluid density [mg/mm^3]
%     FLD.rho = 1e-6;                     % fluid density [mg/mm^3] %MSA: changing temporarily to see the effect of fluid mass

    zz = (-0.5*MP.length_BM:MP.dZ:0.5*MP.length_BM)';
    idx = knnsearch(MP.zz(:),zz(:));

%     FLD.nBot = find(strcmp(Nd.name,MP.BMC)); % BM centerline nodes number
    FLD.nBot = zeros(MP.nz,numel(MP.BMC));
    for i = 1:numel(MP.BMC)
        FLD.nBot(:,i) = find(strcmp(Nd.name,MP.BMC{i})); % BM centerline nodes number
    end
    FLD.nTop = find(strcmp(Nd.name,'BB'));
%     FLD.nTop = [find(strcmp(Nd.name,'E1')),find(strcmp(Nd.name,'CC'))]; %
%     FLD.nTop = reshape(find(strncmp(Nd.name,'E',1)),1201,[]);
%     FLD.nTop = find(strcmp(Nd.name,'BB')); % 
    FLD.nBot = FLD.nBot(idx,:); 
    FLD.nTop = FLD.nTop(idx,:);
    % FLD.beta_wBM = transpose(linspace(0.2,0.8,length(FLD.nBot)));
    FLD.beta_wBM = transpose(linspace(0.2,0.8,length(FLD.nBot)));
    % FLD.beta_wBM = MP.Prop_fit.width_BM(MP.xx)./300;
    FLD.idx = idx;

    global coef
%     coef.dcf_23_A = 0.6;
%     coef.dcf_32_A = 0.6;
%     coef.dcf_23_B = 0.6;
%     coef.dcf_32_B = 0.6;

    % [dcf_32, dcf_23] = dcf_gradient(MP);

    FLD.dcf_32 = 0.6;                   % Dimension conversion factor from 3-D to 2-D (acceleration felt by fluid boudnary)
    FLD.dcf_23 = 0.6;                   % Dimension conversion factor from 2-D to 3-D (pressure felt by CP strutural surfaces)

end

function [dcf_32, dcf_23] = dcf_gradient(MP)

    global coef
    gscale = 'LOG';
    ab_dist = 8000;     % Distance between the reference points at base (z = 2mm) and apex (z = 10 mm)
    if strcmpi(gscale,'LOG')
        grad_32 = log(coef.dcf_32_A./coef.dcf_32_B)/ab_dist;
        grad_23 = log(coef.dcf_23_A./coef.dcf_23_B)/ab_dist;
    else
        grad_32 = (coef.dcf_32_A-coef.dcf_32_B)/ab_dist;
        grad_23 = (coef.dcf_23_A-coef.dcf_23_B)/ab_dist;
    end

    zz = (-0.5*MP.length_BM:MP.dZ:0.5*MP.length_BM)';
    if MP.loc<=6       % values will be referenced to the basal properties at Z = 2mm
        zz = 1e3*(MP.loc-2) + zz;
        dcf_32_0 = coef.dcf_32_B;
        dcf_23_0 = coef.dcf_23_B;
    else                % values will be referenced to the apiaal properties at Z = 10mm
        zz = 1e3*(MP.loc-10) + zz;
        dcf_32_0 = coef.dcf_32_A;
        dcf_23_0 = coef.dcf_23_A;
    end

    if strcmpi(MP.gscale(1:3),'LIN')
        dcf_32 = dcf_32_0 + grad_32*zz;
        dcf_23 = dcf_23_0 + grad_23*zz;
    elseif strcmpi(MP.gscale(1:3),'LOG')
        dcf_32 = dcf_32_0*exp(grad_32*zz);
        dcf_23 = dcf_23_0*exp(grad_23*zz);
    else
        dcf_32 = dcf_32_0 + 0*zz;
        dcf_23 = dcf_23_0 + 0*zz;
    end

end
