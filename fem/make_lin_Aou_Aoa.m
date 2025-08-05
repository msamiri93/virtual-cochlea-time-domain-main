function G = make_lin_Aou_Aoa(G, MP, Nd, El, OHC)

    HB = [OHC.HB];

    ndof = 6;
    HBdir = 1;

    nANK = find(strcmp(El.name,'ANK'));
    nd1 = El.Nd1(nANK); nd2 = El.Nd2(nANK);

    dx1 = (nd1-1)*ndof + HBdir;
    dx2 = (nd2-1)*ndof + HBdir;

    % dxHB = dx2 - dx1
    dpodxx = eval_dpodxx(HB);
    dpodxa = eval_dpodxa(HB);

    odof = MP.dof(5).n;

    row = 1:MP.dof(5).n;

    partsAou = zeros(2*odof,4);
    cntAou = 0;

    [partsAou,cntAou] = fillT(row,dx2,dpodxx,partsAou,cntAou);
    [partsAou,~] = fillT(row,dx1,-dpodxx,partsAou,cntAou);

    adof = MP.dof(4).n;
    col = 1:adof;

    partsAoa = zeros(odof,4);
    cntAoa = 0;

    [partsAoa,~] = fillT(row,col,dpodxa,partsAoa,cntAoa);

    Aou = sparse(partsAou(:,1),partsAou(:,2),partsAou(:,3),odof,Nd.tdof);
    Aoa = sparse(partsAoa(:,1),partsAoa(:,2),partsAoa(:,3),odof,adof);
    % reduce the matrix to unknown degrees of freedom
    G.Aou = -Aou(:,Nd.BC);
    G.Aoa = -Aoa;
    G.Aoo = speye(odof);
end

function dpodx = eval_dpodxx(HB)

    z = [HB.z];
    Xo = [HB.Xo];
    kF = [HB.kF];
    kR = [HB.kR];
    xa = [HB.xa]; % at rest
    xx = [HB.xx]; % at rest
    kBT = 4e-3;
%                 / z (Xo + xa - xx) \
%      kF kR z exp| ---------------- |
%                 \        kBT       /
% ----------------------------------------
%     /            / z (Xo + xa - xx) \ \2
% kBT | kF + kR exp| ---------------- | |
%     \            \        kBT       / /

    dpodx = (kF.*kR.*z.*exp((z.*(Xo + xa - xx))/kBT))./(kBT*(kF + kR.*exp((z.*(Xo + xa - xx))/kBT)).^2);
    
end

function dpodxa = eval_dpodxa(HB)

    z = [HB.z];
    Xo = [HB.Xo];
    kF = [HB.kF];
    kR = [HB.kR];
    xa = [HB.xa]; % at rest
    xx = [HB.xx]; % at rest
    kBT = 4e-3;

%                   / z (Xo + xa - xx) \
%        kF kR z exp| ---------------- |
%                   \        kBT       /
% - ----------------------------------------
%       /            / z (Xo + xa - xx) \ \2
%   kBT | kF + kR exp| ---------------- | |
%       \            \        kBT       / /

    dpodxa = -(kF.*kR.*z.*exp((z.*(Xo + xa - xx))/kBT))./(kBT*(kF + kR.*exp((z.*(Xo + xa - xx))/kBT)).^2);
    
end