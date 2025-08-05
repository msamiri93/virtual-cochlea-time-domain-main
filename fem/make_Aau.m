function Aau = make_Aau(MP,Nd,El,OHC)

    HB = [OHC.HB];
    ndof = 6;   % number of dof for one node
    eANK = find(strcmpi(El.name,'ANK'));
    eANK = eANK([OHC.idx]);
    nd1 = El.Nd1(eANK); nd2 = El.Nd2(eANK);
%     nOHB = find(strcmpi('OHB',El.name));
%     nd1 = El.Nd1(nOHB); nd2 = El.Nd2(nOHB);

    dof = MP.dof(4).n;
    partsA = zeros(2*dof,4); % 2 structural points are used for evaluating adaption at each cross-section
    cntA = 0;
    HBdir = 1;

    % m = diag(reshape([HB.kA].*[HB.kG].*[HB.gamma],dof,1));
    m = diag(reshape([HB.kA].*[HB.kG]*1e-3,dof,1));
    rdof = 1:MP.dof(4).n;
    cdof = (nd1-1)*ndof + HBdir;
    [partsA,cntA] = fillR(rdof,cdof,m,partsA,cntA);
    cdof = (nd2-1)*ndof + HBdir;
    [partsA,cntA] = fillR(rdof,cdof,-m,partsA,cntA);

    Aau = sparse(partsA(:,1),partsA(:,2),partsA(:,3),dof,Nd.tdof);
    % reduce the matrix to unknown degrees of freedom
    Aau = Aau(:,Nd.BC);

end