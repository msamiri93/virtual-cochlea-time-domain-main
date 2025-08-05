function [Aue,cfOHC] = make_Aue(MP,Nd,El,OHC)

    % OHC force implementation, proportional to Vm
    ndof = 6;
    nOHC = find(strcmpi('OHC',El.name));
    nOHC = nOHC([OHC.idx]);
    nd1OHC = El.Nd1(nOHC);
    nd2OHC = El.Nd2(nOHC);

    edof = MP.dof(3).n;

    % circuit model has 5 nodes per section
    cdof_node2 = 2:5:edof; % V2
    cdof_node3 = 3:5:edof; % V3

    % damp the force at the ends
    nz = length(MP.xx);
    nn = 10;
    b = 100;
    dd = dsigmf(MP.xx,[b,MP.xx(nn),b,MP.xx(end-nn)]);

    % motility inactivation
    
    global batch_coef
    if isfield(batch_coef,'iloc')
        iloc = knnsearch(MP.xx(:),batch_coef.iloc); % in units of [10um]
        ispan = 40; % in units of [10um]

        ds = zeros(size(dd));
        xx = MP.xx(iloc-0.5*ispan):1e-3*MP.dZ:MP.xx(iloc+0.5*ispan);
        ds((iloc-0.5*ispan):(iloc+0.5*ispan)) = 0.5*(1+cos(2*pi*(xx-MP.xx(iloc))/(ispan*1e-3*MP.dZ)));
        
        % ds = dsigmf(MP.xx,[b,MP.xx(iloc-0.5*ispan),b,MP.xx(iloc+0.5*ispan)]);
        
        dd = dd - ds;
    end

    M = [OHC.M];
    % cfOHC = reshape([OHC.count].*[M.gain_fOHC],nz,1); % mind the sign
    cfOHC = reshape([OHC.count].*[M.gain_fOHC].*dd,nz,1); % mind the sign
    fdir = El.dir(nOHC,:);

    num_entries = 2*size(fdir,2)*2*length(cdof_node2); % number of values going into sparse matrix, calculated based on how many entries are evaluated
    partsA = zeros(num_entries,4);
    cntA = 0;

    for idd = 1:size(fdir,2)
        m = diag(fdir(:,idd).*cfOHC);

        rdof = (nd1OHC-1)*ndof + idd;      
        [partsA,cntA] = fillR(rdof,cdof_node2,-m,partsA,cntA);
        [partsA,cntA] = fillR(rdof,cdof_node3,m,partsA,cntA);

        rdof = (nd2OHC-1)*ndof + idd;
        [partsA,cntA] = fillR(rdof,cdof_node2,m,partsA,cntA);
        [partsA,cntA] = fillR(rdof,cdof_node3,-m,partsA,cntA);
    end

    Aue = sparse(partsA(:,1),partsA(:,2),partsA(:,3),Nd.tdof,edof);
    % reduce the matrix to unknown degrees of freedom
    Aue = Aue(Nd.BC,:);

end