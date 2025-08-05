function FLD = set_FEfluid_BCs(MP,FLD)

    p = FLD.Nd.x;
    h0 = FLD.h0;
    HOC = FLD.HOC;
    L = FLD.L;
    nnode = length(p);
    FLD.name = cell(nnode,1);
    FLD.name(:) = {''};      

    indTM = FLD.Nd.bnd.TOP_SURFACE;
    indBM = FLD.Nd.bnd.BOTTOM_SURFACE;
    
    % indTM = find(p(:,1) < L + h0/100 & p(:,1) > - h0/100 & p(:,2) > 0 & p(:,2) < HOC + h0/100);
    % indBM = find(p(:,1) < L + h0/100 & p(:,1) > - h0/100 & p(:,2) < 0 & p(:,2) > -HOC - h0/100);
    
    FLD.name(indTM,1) = {'TM'};
    FLD.name(indBM,1) = {'BM'}; 
    
    % indOW = find(abs(p(:,1)) < h0/100 & p(:,2) > 0 );
    % indRW = find(abs(p(:,1)) < h0/100 & p(:,2) < 0 );

    indOW = FLD.Nd.bnd.OW_SURFACE;
    indRW = FLD.Nd.bnd.RW_SURFACE;
    
    % for mesh mesh_12mm_cochlea_quad_5.mat % extracted once to be reused.
%     indOW = [3,5026,5027,5028,5029,5030,5031,5032,5033,5034,5035,5036,5037,...
%         5038,5039,5040,5041,5042,5043];
%     indRW = [4,5065,5066,5067,5068,5069,5070,5071,5072,5073,5074,5075,5076,...
%         5077,5078,5079,5080,5081,5082];

    FLD.name(indOW,1) = {'OW'};
    FLD.name(indRW,1) = {'RW'};  
    
    if MP.opt_wtr_shnt == 1
        ind_TS = find(p(:,1) > 0.95*L + h0/100 & p(:,1) < L + h0/100 & p(:,2) > 0 & p(:,2) < HOC + h0/100);
        ind_BS = find(p(:,1) > 0.95*L + h0/100 & p(:,1) < L + h0/100 & p(:,2) < 0 & p(:,2) > -HOC - h0/100);
        FLD.ind_TS = ind_TS;
        FLD.ind_BS = ind_BS;
    end

    FLD.indTM = indTM;
    FLD.indBM = indBM;

    FLD.indOW = indOW;
    FLD.indRW = indRW;

    FLD.BC = true(nnode,1);
    FLD.BC(strncmp(FLD.name,'OW',2) | strncmp(FLD.name,'RW',2)) = false;
%     FLD.BC(strncmp(FLD.name,'OW',2)) = false;
    FLD.tdof = nnode;
    FLD.rdof = sum(FLD.BC);
    
end