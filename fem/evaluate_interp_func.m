function El = evaluate_interp_func(Nd,El,xi,eta)

%     global GaussQuad;    
%     xi = GaussQuad.xi;
%     eta = GaussQuad.eta;    
    
    NQ = length(xi);
    type = El.type;

    elm_psi = zeros(type,NQ,El.N);
    elm_gpsi = zeros(2*type,NQ,El.N);
    elm_hs = zeros(1,NQ,El.N);
    
    for ie = 1:El.N       
        enodes = El.node(ie,:);
        xx =  Nd.x(enodes,:);        
        for iq = 1:NQ            
            [psi, gpsi, hs] = elm_interp(type, xx, xi(iq), eta(iq));            
            elm_psi(:,iq ,ie) = psi';
            elm_gpsi(:,iq,ie) = gpsi(:)';
            elm_hs(:,iq,ie) = hs;            
        end           
    end
    El.elm_psi = elm_psi;
    El.elm_gpsi = elm_gpsi;
    El.elm_hs = elm_hs;    
end