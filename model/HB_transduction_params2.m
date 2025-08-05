function HB = HB_transduction_params2(MP, dZ)
% %      First version written: 2009.10.13
% %      Last Modified: 2013.11.09
% %      by Jong-Hoon Nam
% %      at University of Rochester   
% % -----------------------------------------------------------------------
%
% Define OHC transduction channel parameters
% loc = location a: apex, b: base, number: 0(basal end)-1(apical end)
%
% Units are um, pN, ms, aJ
% Ks = 1 mN/m = 1000 pN/um, stereocilia stiffness
% c = 100 nN-s/m = 100 pN-ms/um, damping of HB
% m = 150 ng, mass carried by a single HB
% Xo = 0 nm = 0.000 um, set point to define resting open probability
% Kg = 4 mN/m = 4000 pN/um, Kg = N*gamma^2*(stiff of gating spring)
% N = 60, number of MT channels per stereocilia
% b = 2 nm = 0.002 um, gating swing
% gamma = 0.11, geometric gain of HB (=dxTL/dxHB)
% z = 1.2 pN, z = gamma*(Kg/(N*gamma^2))*b
% XCa = 50 nm = 0.05 um, shift of po-X curve due to Ca-bind
% C0 = 0.1 uM, [Ca2+] at the channel site when channel remained closed
% C1 = 50 uM, [Ca2+] at the channel site when channel remained open
% KD0 = 5 uM, Ca dissoc. coefficient from channel when channel is closed
% KD1 = 50 uM, Ca dissoc. coefficient from channel when channel is open
% kB = 0.5 ms^-1uM^-1, Ca binding coeff
% kF = 10 ms^-1, channel forward rate coeff.
% kR = 10 ms^-1, channel backward rate coeff.
% kA = 0.25 ms^-1, slow adaptation rate
%
% Version 2 receives generalized longitudinal locations as an input.

    loc = MP.loc;
    A = set_HB_props_at_two_locations(MP,'a');%%%%%%%%
    B = set_HB_props_at_two_locations(MP,'b');%%%%%%%%%%%

    if ~exist('dZ','var')
        dZ = 0;
    end

    HB = set_HB_props_at_arb_locations(loc, dZ, A, B);
    
    HB.cState = 0.1*ones(2,5);      % |C0, C1, C1, C2, C4|
                                    % |O0, O1, O2, O3, O4|  
                                    
    HB.Fo = 0;
    HB.XA = 0e-3;      
    
end

function HB = set_HB_props_at_two_locations(MP,loc)%%%%%%%%%%%%%

    if strcmpi(loc(1),'A')
        HB.Ks = 5e3;      % stiffness of single HB
        HB.c = 200;       % should be ~50, based on TM width 210 um half of which is the frictional surface.
                          % In the whole CP model, this value is neglected and replaced with a direct calculation in make_C0. 
                          % In make_CO calculated from the sub-TM viscous friction.
        
        HB.m = 0.001;
                        
        HB.Xo = 4*1.5e-3;        % this determines the initial state, increasing Xo decreases po,rest 
                                % Xo = -HB.N*HB.z*po,rest/k_OHB

        HB.N = 60;          % number of tip links in HB        
        HB.b = 1.54*1e-3;      % gaing swing
        % HB.b = coef.HB_b_A*0.8e-3; % [um]

        HB.gamma = 0.1;        
        
        HB.kG = 1.0*6e3;
        HB.KG = HB.N*HB.gamma^2*HB.kG;       % Kg = N*gamma^2*(stiff of actual TLA)        
        HB.z = HB.gamma*HB.b*HB.kG;
        % HB.kES = 5*0.6e3;
        % HB.kES = 0.14*HB.gamma*HB.kG;
        HB.kES = 0.1*HB.gamma*HB.kG;
        HB.k = HB.z/8;
        % HB.po = 0.36;
        HB.XCa = 3.0e-3;    % greater XCa results in stronger effect of Ca2+ on fast adaptation

        HB.C0 = 1;
        HB.C1 = 40;

        HB.KD0 = 40;
        HB.KD1 = 1;

        HB.kB = 0.4;
               
        HB.kF = 6*10;
        HB.kR = 6*10;
        
        % HB.kA = 6e-3;
        % HB.kA = 0.05e-3; % Prodanovic
        % HB.kA = 0.1;
        HB.kA = 5e-2;

        HB.sigs = 0.5*[HB.KD0*HB.C0, HB.KD1*HB.C1, HB.kF];
        
    elseif strcmpi(loc(1),'B')
        
        HB.Ks = 40e3;           % This value is not used        
        HB.c = 100;             % This value is not used
        HB.m = 3.3;             % 3.33 ng = (10 um x 20 um x 50 um x 1e-3 ng/um^3)/(3 HBs)
        HB.m = 0.001;           % Is this used?    
        HB.Xo = 7*0.4e-3;        % this determines the initial state, increasing Xo decreases po,rest  
                
        HB.N = 75;          % number of tip links in HB        
        
        HB.b = 1.2*1e-3;
        % HB.b = coef.HB_b_B*1.0e-3;
               
        % HB.gamma = 0.25;         % geometric gain   
        HB.gamma = 0.25;
        % HB.kG = 6e3;
        HB.kG = 1.0*6e3;
        HB.KG = HB.N*HB.gamma^2*HB.kG;       % Kg = N*gamma^2*(stiff of actual TLA)
        % HB.kES = 5*1.5e3;
        % HB.kES = 0.05*HB.gamma*HB.kG;
        HB.kES = 0.1*HB.gamma*HB.kG;
        HB.z = HB.gamma*HB.b*HB.kG;
        HB.k = HB.z/8;
%         HB.po = 0.36;
        HB.XCa = 0.7e-3;   % greater XCa results in stronger effect of Ca2+ on fast adaptation fCa = gamma*kG*b
        
        HB.C0 = 1.0;
        HB.C1 = 100;
        HB.KD0 = 100;
        HB.KD1 = 1.0;
               
        HB.kB = 0.4;    % default 0.2
      
        HB.kF = 100;
        HB.kR = 100;
        
        % HB.kA = 5*2e-3;
        % HB.kA = 0.05e-3; % Prodanovic
        % HB.kA = 0.05;
        HB.kA = 5e-2;      
                
        HB.sigs = 0.5*[HB.KD0*HB.C0, HB.KD1*HB.C1, HB.kF];        
        
    else
        fprintf(1,'\nUndefined option\n');
        
    end
    
    %     HB.ch_noise = rand(1,3)-0.5;
    
end


function HB = set_HB_props_at_arb_locations(loc, dZ, A, B)
    
    if loc <= 6
        HB0 = B;
    elseif loc > 6
        HB0 = A;
    end
    
       
    PROPS = {'Ks', 'c', 'm', 'Xo', 'N', 'b', 'gamma', 'kG', 'KG', 'z', 'XCa', ...
                'C0', 'C1', 'KD0', 'KD1','kB', 'kF', 'kR', 'kA','sigs','k','kES'};
    
    for ii=1:length(PROPS)
                
        prop = PROPS{ii};
        if isfield(HB0,prop)
            
% %             grad = set_longi_gradient(A,B,prop,'LIN');
% %             HB.(prop) = HB0.(prop)+grad*dZ;        
            
            grad = set_longi_gradient(A,B,prop,'LOG');
            HB.(prop) = HB0.(prop).*exp(grad.*dZ);
        end
    end
    
end