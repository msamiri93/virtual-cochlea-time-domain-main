function FLD = set_FEfluid_BCs_MC(MP,FLD)

    p = FLD.Nd.x;
    h0 = FLD.h0;
    Xo = FLD.Xo;
    Lo = FLD.Lo;
    Ho = FLD.Ho;
    Xs = FLD.Xs;
    Ls = FLD.Ls;
    Hs = FLD.Hs;
    H = FLD.H;
    HOC = FLD.HOC;
    L = FLD.L;
    nnode = length(p);
    FLD.name = cell(nnode,1);
    FLD.name(:) = {''}; 

%     TC_indW = find(abs(p(:,1)+1000)<h0/10 & p(:,2)>-1);
%     TC_indE = find(abs(p(:,1)-2000)<h0/10 & p(:,2)>-1);
%     TC_indS = find(abs(p(:,2))<h0/10);
%     TC_indN = find(abs(p(:,2)-1000)<h0/10);
%     BC_indW = find(abs(p(:,1)+2000)<h0/10 & p(:,2)<-100);
%     BC_indE = find(abs(p(:,1)-3000)<h0/10 & p(:,2)<-100);
%     BC_indS = find(abs(p(:,2)+1200)<h0/10);
%     BC_indN = find(abs(p(:,2)+200)<h0/10);
%     SL_indW = find(abs(p(:,1))<h0/10 & p(:,2)<h0/10 & p(:,2)>h0/10-200);
%     SL_indE = find(abs(p(:,1)-1000)<h0/10 & p(:,2)<h0/10 & p(:,2)>h0/10-200);
%     SL_indS = find(abs(p(:,2)+200)<h0/10 & p(:,1)>h0/10 & p(:,1)<h0/10+1000);
%     SL_indN = find(abs(p(:,2))<h0/10 & p(:,1)>h0/10 & p(:,1)<h0/10+1000);

    indTM = find(p(:,1) < L + h0/100 & p(:,1) > - h0/100 & p(:,2) > 0 & p(:,2) < HOC + h0/100);
    indBM = find(p(:,1) < L + h0/100 & p(:,1) > - h0/100 & p(:,2) < 0 & p(:,2) > -HOC - h0/100);
    
    FLD.name(indTM,1) = {'TM'};
    FLD.name(indBM,1) = {'BM'};

    TC_indW = find(abs(p(:,1)-(Xo-1000))<h0/10 & p(:,2)>H/2+Ho/2);
    TC_indE = find(abs(p(:,1)-(Xo+Lo+1000))<h0/10 & p(:,2)>H/2+Ho/2);
%     TC_indS = find(abs(p(:,2)-(H/2+Hs))<h0/10);
%     TC_indN = find(abs(p(:,2)-(H/2+Hs+MP.topH))<h0/10);
    BC_indW = find(abs(p(:,1)-(Xs-1000))<h0/10 & p(:,2)<-H/2-Hs);
    BC_indE = find(abs(p(:,1)-(Xs+Ls+1000))<h0/10 & p(:,2)<-H/2-Hs);
    BC_indS = find(abs(p(:,2)-min(p(:,2)))<h0/10);

    FLD.name(TC_indW,1) = {'TCW'};
    FLD.name(TC_indE,1) = {'TCE'};    
    FLD.name(BC_indW,1) = {'BCW'};
    FLD.name(BC_indE,1) = {'BCE'};    
    FLD.name(BC_indS,1) = {'BCS'};

    FLD.indTM = indTM;
    FLD.indBM = indBM;
    FLD.indTCW = TC_indW;
    FLD.indTCE = TC_indE;    
    FLD.indTBW = BC_indW;
    FLD.indTBE = BC_indE;  
    FLD.indTBS = BC_indS;  


    FLD.BC = true(nnode,1);
    FLD.BC(strncmp(FLD.name,'TCW',3) | strncmp(FLD.name,'TCE',3) | ...
        strncmp(FLD.name,'BCS',3) | strncmp(FLD.name,'BCW',3)) = false;
    FLD.tdof = nnode;
    FLD.rdof = sum(FLD.BC);

end