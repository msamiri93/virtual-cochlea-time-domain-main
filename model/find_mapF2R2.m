function Nd = find_mapF2R2(Nd)
% %       function [mapF2R,mapR2F,rdof]=find_mapF2R(Nd)
% % -----------------------------------------------------------------------
% % 
% %      Description: generates vectors related index of full dof and reduced dof
% %                       
% %      First written: 2003.07.18
% %      Last modified: 2007.03.27
% %      by Jong-Hoon Nam
% % 

    NDOF = 6;
    tdof = Nd.N*NDOF;
    
    mapR2F = (1:tdof)';

%     mapR2F(1:3) = 0;                            % fix node 1
%     mapR2F(NDOF*(7-1)+1:NDOF*(7-1)+3) = 0;      % fix node 7    
%     mapR2F(3:NDOF:tdof) = 0;                     % fix z-dof
%     
%     mapR2F(4:NDOF:tdof) = 0;
%     mapR2F(5:NDOF:tdof) = 0;
% %     mapR2F(6:NDOF:tdof) = 0;
%     mapR2F(tdof) = 0;

    for ind = 1:Nd.N,
        jdof = (ind-1)*NDOF;
        mapR2F(jdof+1:jdof+NDOF) = Nd.BC(ind,:)'.*(jdof+1:jdof+NDOF)';        
    end
    
    mapR2F(mapR2F==0)=[];
    

    rdof=length(mapR2F);    
    mapF2R=zeros(tdof,1);
    
    for ir=1:rdof,
        mapF2R(mapR2F(ir))=ir;
    end
    
    Nd.mapF2R = mapF2R;
    Nd.mapR2F = mapR2F;
    Nd.tdof = tdof;
    Nd.rdof = rdof;
end