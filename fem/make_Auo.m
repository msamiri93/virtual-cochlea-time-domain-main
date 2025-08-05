function [Auo,cfMET] = make_Auo(MP,Nd,El,OHC,HB)

    % MET force implementation, proportional to open probability
    ndof = 6;
    nOHB = find(strcmpi('OHB',El.name));
    nOHB = nOHB([OHC.idx]);
    nANK = find(strcmpi('ANK',El.name));
    nANK = nANK([OHC.idx]);
    
    nd1OHB = El.Nd1(nOHB);
    nd2OHB = El.Nd2(nOHB);

    odof = MP.dof(5).n;

    col = 1:MP.dof(5).n;

    % damp the force at the ends
    nn = 10;
    b = 100;
    dd = dsigmf(MP.xx,[b,MP.xx(nn),b,MP.xx(end-nn)]);

    cfMET = reshape([OHC.count].*[HB.N].*[HB.z],odof,1);
    fdir = El.dir(nANK,:);

    partsA = zeros(2*size(fdir,2)*odof,4);
    cntA = 0;

%     ------> + dx | + p_o
%        _
%       | | -----> f_MET = nd2OHB
%       | |
%       | |
%       | |
%       | |
%       | |
%       \ / <----- f_MET = - nd1OHB
%
    for idd = 1:size(fdir,2)
        m = diag(fdir(:,idd).*cfMET);

        row = (nd1OHB-1)*ndof + idd;      
        [partsA,cntA] = fillR(row,col,m,partsA,cntA);

        row = (nd2OHB-1)*ndof + idd;
        [partsA,cntA] = fillR(row,col,-m,partsA,cntA);
    end

    Auo = sparse(partsA(:,1),partsA(:,2),partsA(:,3),Nd.tdof,odof);
    % reduce the matrix to unknown degrees of freedom
    Auo = Auo(Nd.BC,:);

end