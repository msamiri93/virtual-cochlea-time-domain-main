function [S, M] = OHC_electrical_params2(MP,dZ)

    % S: stereociliar bundle
    % M: cell body
    
   if ~exist('dZ','var')      
        [S, M] = set_OHC_electrical_props_at_arb_locations(MP, 0);
    else
        [S, M] = set_OHC_electrical_props_at_arb_locations(MP, dZ);
    end
    
end


function [S, M] = set_OHC_electrical_props_at_arb_locations(MP, dZ)

    loc = MP.loc;
    [Sa, Ma, Sb, Mb] = set_OHC_electrical_props_at_two_locations(MP);
    
    if loc <= 6
        M0 = Mb;
        S0 = Sb;
    else
        M0 = Ma;
        S0 = Sa;
    end
    
    

    prop = 'Gmax';
    grad = set_longi_gradient(Ma,Mb,prop);
    M.(prop) = M0.(prop)*exp(grad*dZ);
    
    prop = 'VG05';
    grad = set_longi_gradient(Ma,Mb,prop);
    M.(prop) = M0.(prop)*exp(grad*dZ);
    
    prop = 'eK';
    grad = set_longi_gradient(Ma,Mb,prop);
    M.(prop) = M0.(prop)*exp(grad*dZ);
    
    prop = 'Qmax';
    grad = set_longi_gradient(Ma,Mb,prop);
    M.(prop) = M0.(prop)*exp(grad*dZ);
    
    prop = 'CLin';
    grad = set_longi_gradient(Ma,Mb,prop);
    M.(prop) = M0.(prop)*exp(grad*dZ);
    
    prop = 'VC05';
    grad = set_longi_gradient(Ma,Mb,prop);
    M.(prop) = M0.(prop)*exp(grad*dZ);
    
    prop = 'zC';
    grad = set_longi_gradient(Ma,Mb,prop);
    M.(prop) = M0.(prop)*exp(grad*dZ);
    
    prop = 'tauK';
    grad = set_longi_gradient(Ma,Mb,prop);
    M.(prop) = M0.(prop)*exp(grad*dZ);
    
    prop = 'gain_fOHC';
    grad = set_longi_gradient(Ma,Mb,prop);
    M.(prop) = M0.(prop)*exp(grad*dZ);

    prop = 'Ek';
    grad = set_longi_gradient(Ma,Mb,prop);
    M.(prop) = M0.(prop)*exp(grad*dZ);
    
    prop = 'V05';
    grad = set_longi_gradient(Ma,Mb,prop);
    M.(prop) = M0.(prop)*exp(grad*dZ);
        
    prop = 'Vs';
    grad = set_longi_gradient(Ma,Mb,prop);
    M.(prop) = M0.(prop)*exp(grad*dZ);
    
    prop = 'Gmax';
    grad = set_longi_gradient(Sa,Sb,prop);
    S.(prop) = S0.(prop)*exp(grad*dZ);
    
    prop = 'C';
    grad = set_longi_gradient(Sa,Sb,prop);
    S.(prop) = S0.(prop)*exp(grad*dZ);

end

function [Sa, Ma, Sb, Mb] =  set_OHC_electrical_props_at_two_locations(MP)

        if MP.vMC == 1
            Gm = 1e4;
            Gs = 1e4;
            Ek = 0;
        else
            Gm = 1;
            Gs = 1;
            Ek = 1;
        end

%        case 'a'                % at 0.7 kHz
%             fprintf(1,['\nOHC prob set [' loc '] called...\n']);
            % M.Gmax = Gm*45;        % max cell membrane conductance in nS
			M.Gmax = Gm*39;
            M.VG05 = -70;       % half maximal MP for membrane conductance in mV
            M.eK = 18;          % slope width of M.G in mV
            M.Qmax = 1.5;       % maximum nonlinear mobile charge in pC
            M.CLin = 15;        % linear component of cell membrane capacitance in pF
            M.VC05 = -40;       % MP at maximal capacitance in mV
            M.zC = 0.97;        % valence of mobile membrane charge
            M.tauK = 5;         % time const of membrane potassium channel in ms
            
            M.gain_fOHC = 100;  % gain for fOHC in pN/mV 

            M.Ek = Ek*75;          %OHC equilibrium potential in mV, source: Mistrik 2009
            
            M.V05 = -40;        % half maximal membrane voltage in mv
            M.Vs = 28;          % slope of fOHC-mV curve in mV
            % fOHC = gain_fOHC*Vs*(-1/2 + 1/(1+exp(-(Vm-V05)/Vs)))
            
            S.Gmax = Gs*27;        % stereociliar conductance in nS
            % S.C = 1/6*M.CLin;   % stereociliar capacitance in pF
			S.C = 9;
                        
        Ma = M;
        Sa = S;
            
%        case 'b'                % about 20 kHz
%             fprintf(1,['\nOHC prob set [' loc '] called...\n']);

            % M.Gmax = Gm*320;       % max cell membrane conductance in nS 
            M.Gmax = Gm*320;
            M.VG05 = -70;       % half maximal MP for membrane conductance in mV
            M.eK = 18;          % slope width of M.G in mV
            M.Qmax = 0.6;       % maximum nonlinear mobile charge in pC
            M.CLin = 4.3;       % linear component of cell membrane capacitance in pF
            
            M.VC05 = -40;       % MP at maximal capacitance in mV
            M.zC = 0.94;        % valence of mobile membrane charge
            M.tauK = 0.5;       % time const of membrane potassium channel in ms
            
            M.gain_fOHC = 100;  % gain for fOHC in pN/mV

            M.Ek = Ek*75;          %OHC equilibrium potential in mV, source: Mistrik 2009

            M.V05 = -40;        % half maximal membrane voltage in mv
            M.Vs = 28;          % slope of fOHC-mV curve in mV
            % fOHC = gain_fOHC*Vs*(-1/2 + 1/(1+exp(-(Vm-V05)/Vs)))           
            
            S.Gmax = 90;       % stereociliar conductance in nS
            % S.C = 1/6*M.CLin;   % stereociliar capacitance in pF 
			S.C = 2.6;
                        
% %             M.Gmax = 300;       % max cell membrane conductance in nS                        
% %             S.Gmax = 90;       % stereociliar conductance in nS

            S.Gmax = Gs*90;       % stereociliar conductance in nS
            
        Mb = M;
        Sb = S;
 
    
end