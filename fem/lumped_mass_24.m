function  El = lumped_mass_24(Nd, El,  MP)

    eBMa = strcmp('BMa',El.name);
    eBMp = strcmp('BMp',El.name);
    eBM = eBMp;

    eBM = find(eBM);

    thick_eX = thickness_BM(Nd, El, eBM, MP);
   
    eX = eBM; mthick = thick_eX;        
    mfactor = 1.5*mthick./El.dim(eX,2); % factor 1.5 in order to include supporting cells in OC                  
    El.rho(eX) = mfactor.*El.rho(eX);        

end

function thick_eX = thickness_BM(Nd, El, eX, MP)

    Prop = MP.Prop_fit;
    if mean(El.type(eX))>2.5,   % if it is a plate
        nds = [El.Nd1(eX), El.Nd2(eX), El.Nd3(eX), El.Nd4(eX)];
    else                        % else if it is a beam
        nds = [El.Nd1(eX), El.Nd2(eX)];
    end
    xi = mean(Nd.X(nds),2);
    zi = mean(Nd.Z(nds),2);

    a = 1e-3; % convert [um] to [mm] for interpolation
    peak_thick = feval(Prop.('thick_BMz'),zi*a); % depth  

    thick_fluid = feval(Prop.('thick_fluid'),zi*a);
    % fluid-mass gradient may be ignored (at the singe stimulating frequency, the fluid mass may not increase along the length)
    
    e_BMpwidth = feval(Prop.('width_BMP'),zi*a);
    e_BMawidth = feval(Prop.('width_BMA'),zi*a);
    xi = xi./(e_BMpwidth + e_BMawidth); 
    xap =  e_BMawidth./(e_BMpwidth + e_BMawidth);

    x_profile = zeros(size(xi));
    za = contains(El.name(eX),'a'); % zona arcuate
    zp = contains(El.name(eX),'p'); % zona pectinate
    ze = xi==0 | xi==1; % edge
    % x_profile(za) = -36*xi(za).*(xi(za)-xap(za)); 
    % x_profile(zp) = -9*(xi(zp)-xap(zp)).*(xi(zp)-1);
    x_profile(za) = -25*xi(za).*(xi(za)-xap(za));
    x_profile(zp) = -6.25*((xi(zp)-xap(zp)).^0.6).*(xi(zp)-1);   
    thick_eX = x_profile.*peak_thick;
    % thick_eX(ze) = 0.1;

    thick_eX = thick_eX + 0.5*thick_fluid;
end