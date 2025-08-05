function El = model_moduli5_24a(El, Nd, MP, opt_lgradient)
    
    % this can also be included in the mat file instead of hard-coding
    eBMa = strcmp('BMa',El.name);   % radial BM fiber layer arcuate zone
    eBMp = strcmp('BMp',El.name);   % radial BM fiber layer BM pectinate zone
    eBMz = strcmp('BMza',El.name) | strcmp('BMzp',El.name);   %   
    
    eOHC = strcmp('OHC',El.name);
    eIPC  = strncmp('IPC',El.name,3);
    eOPC  = strncmp('OPC',El.name,3);
    eDCb = strncmp('DCb',El.name,3); % changed to account for Dbc
    eDCp = strcmp('DCp',El.name);
    eDCz = strcmp('DCz',El.name);
    eDCr = strcmp('DCr',El.name);
    eRLx1 = strcmp('RLx1',El.name);
    eRLz = strcmp('RLz',El.name);
    eRLp = strcmp('RLp',El.name);
    eRLpz = strcmp('RLpz',El.name);
    eTMx1 = strcmp('TMx1',El.name); % TM attachment
    eTMx2 = strcmp('TMx2',El.name); % TM body
    eTMz0 = strcmp('TMz0',El.name);
    eTMz1 = strcmp('TMz1',El.name);
    eTMz2 = strcmp('TMz2',El.name);
    eOHB = strcmp('OHB',El.name); % OHC bundle and linear spring Nd1:bb, Nd2:e1
    exHB = strcmp('xHB',El.name); % Pseudo rod at OHC bundle location
    eANK = strcmp('ANK',El.name); % Pseudo link to conveniantly calculate translational displ.
    eTCz = strcmp('TCz',El.name);

    El.type = ones(El.N,1)*0.1; % Default Timoshenko beam: 0; Euler beam: 0.1
    El.type(eRLz) = 0.1;
    El.type(eDCz) = 2;
    El.type(eANK) = 2;        
    El.type(eTMx2) = 0.1;    
    El.type(eTCz | eRLz | eRLpz) = 2;
    % element type 0.2 means tapered beam, -:Nd1, +:Nd2, to define the thin end
%     El.type(eDCb) = -0.2;

    El.oL = zeros(El.N,1);
    El.L = zeros(El.N,1);
    El.odir = zeros(El.N,3);
    El.dir = zeros(El.N,3);
    
    El = get_vector(El, Nd, 0);
    
    El.dim = zeros(El.N,2);  % dim = [bx, bz], note that x is the primary axis    
    El = set_dimension_gradient(El, Nd, MP, opt_lgradient);    
    
    % Elements with rectangular cross-section
    eRectangular = eBMa | eBMp | eBMz | eRLx1 | eRLz | eRLp | eRLpz | eTMx1 | eTMx2 | eTMz0 | eTMz1 | eTMz2 | eOHB | exHB | eANK;
    % Elements with circular cross-section
    eCircular = eOHC | eIPC | eOPC | eDCb | eDCp | eDCz | eDCr | eTCz;
    
    El.A = zeros(El.N,1);
    El.Iz = zeros(El.N,1);
    El.Iy = zeros(El.N,1);
    
    El.A(eRectangular) = El.dim(eRectangular,1).*El.dim(eRectangular,2);
    El.A(eCircular) = (pi/4)*El.dim(eCircular,1).^2;
    
    El.Iz(eRectangular) = (1/12)*El.dim(eRectangular,1).*El.dim(eRectangular,2).^3;
    El.Iz(eCircular) = (pi/64)*El.dim(eCircular,1).^4;
    
    El.Iy(eRectangular) = (1/12)*El.dim(eRectangular,2).*El.dim(eRectangular,1).^3;
    El.Iy(eCircular) = (pi/64)*El.dim(eCircular,1).^4;

    El.A(eOHC | eDCb) = 3*El.A(eOHC | eDCb);
    El.Iz(eOHC | eDCb) = 3*El.Iz(eOHC | eDCb) + 2*El.A(eOHC | eDCb).*El.dim(eOHC | eDCb,1);
    El.Iy(eOHC | eDCb) = 3*El.Iz(eOHC | eDCb);

    % Assigned these values for every element for initialization, if any
    % element has these values, it means that it was not assigned a proper MP.
    El.YM = 123e6*ones(El.N,1);  
    El.SM = 12.3e6*ones(El.N,1);
    El.nu = 0.3*ones(El.N,1);

    Prop = MP.Prop_fit;
    P_names = fieldnames(Prop);
    idx = find(strncmp(P_names,'Y_',2));

    YM_names = replace(P_names(idx),'Y_','');

    a = 1e-3; % 1e-3 coef is to translate [um] to [mm] for interpolation

    for i = 1:numel(YM_names)
        eX = contains(El.name,YM_names{i});
        zi = mean([Nd.Z(El.Nd1(eX)),Nd.Z(El.Nd2(eX))], 2);
        El.YM(eX) = feval(Prop.(P_names{idx(i)}),zi*a);
        El.SM(eX) = 1/3*El.YM(eX);
    end

    % hinge ratio condition
    eX = contains(El.name,'DCb');
    zi = mean([Nd.Z(El.Nd1(eX)),Nd.Z(El.Nd2(eX))], 2);
    El.YM(eX) = feval(Prop.('DC_hinge'),zi*a).*El.YM(eX);
    El.SM(eX) = 1/3*El.YM(eX); 

    eX = contains(El.name,'DCbc');
    zi = mean([Nd.Z(El.Nd1(eX)),Nd.Z(El.Nd2(eX))], 2);
    El.YM(eX) = (1-feval(Prop.('DC_hinge'),zi*a)).*El.YM(eX);
    El.SM(eX) = 1/3*El.YM(eX);

    El.rho = 1.0e-3*ones(El.N,1); % nanogram/um^3

    El = lumped_mass_24(Nd, El,  MP);

%     El.rho(eTMx1 | eTMx2 | eTMz0 | eTMz1 | eTMz2) = 0.5e-3;
    El.rho(eTMx1 | eTMx2) = 1e-6;
    El.rho(eTMz0 | eTMz1) = 1e-6;

    eZ = contains(El.name,'z') & ~contains(El.name,'BM');
    El.rho(eZ) = 1e-6;

    % eZ = contains(El.name,'z');
    % El.rho(eZ) = 1e-6;

    NDOF = 6;
    El.fi = zeros(El.N,2*NDOF);    

end

function El = set_dimension_gradient(El, Nd, MP, opt_lgradient)

    Prop = MP.Prop_fit;

    if isfield(MP,'nz')
        nz = MP.nz;
    else
        nz = length(MP.xx);
    end

    if opt_lgradient
        xx = MP.xx(:); % in mm
    else
        xx = MP.loc*ones(nz,1);
    end

    eBMa = strcmp('BMa',El.name);   % BM arcuate zone
    eBMp = strcmp('BMp',El.name);   % BM pectinate zone
    % eBMz = strcmp('BMza',El.name) | strcmp('BMzp',El.name);
    
    eOHC = strcmp('OHC',El.name);
    eIPC  = strncmp('IPC',El.name,3);
    eOPC  = strncmp('OPC',El.name,3);
    eDCb = strncmp('DCb',El.name,3);
    eDCp = strcmp('DCp',El.name);
    eDCz = strcmp('DCz',El.name);
    eDCr = strcmp('DCr',El.name);
    eRLx1 = strcmp('RLx1',El.name);
    eRLp = strcmp('RLp',El.name);
    eRLz = strcmp('RLz',El.name);
    eRLpz = strcmp('RLpz',El.name);
    eTMx1 = strcmp('TMx1',El.name);
    eTMx2 = strcmp('TMx2',El.name);
    eTMz0 = strcmp('TMz0',El.name);
    eTMz1 = strcmp('TMz1',El.name);  
%     eTMz2 = strcmp('TMz2',El.name);
    eOHB = strcmp('OHB',El.name);
    exHB = strcmp('xHB',El.name);
    eANK = strcmp('ANK',El.name);
    eTCz = strcmp('TCz',El.name);
            
    eBMza = strcmp('BMza',El.name);     % [YJ updated,11/01/13]
    eBMzp = strcmp('BMzp',El.name);     % [YJ updated,11/01/13]
    
    El.dim = zeros(El.N,2);  % dim = [by, bz], note that x is the primary axis
                             %      bz is the thickness in primary bending
                             %      direction
    
    a = 1e-3; % 1e-3 coef is to translate [um] to [mm] for interpolation
    eX = eBMa; bz = 'thick_BMp';
    zi = mean([Nd.Z(El.Nd1(eX)),Nd.Z(El.Nd2(eX))], 2);
    El.dim(eX,1) = MP.dZ;
    El.dim(eX,2) = 0.333*feval(Prop.(bz),zi*a); % half BMp thickness, can be applied separately in GP explicitly if needed

    eX = eBMp;  bz = 'thick_BMp';
    zi = mean([Nd.Z(El.Nd1(eX)),Nd.Z(El.Nd2(eX))], 2);
    El.dim(eX,1) = MP.dZ;
    El.dim(eX,2) = feval(Prop.(bz),zi*a);    

    eX = eBMzp; bz = 'thick_BMp';       
    zi = mean([Nd.Z(El.Nd1(eX)),Nd.Z(El.Nd2(eX))], 2);
    El.dim(eX,1) = 0.333*feval(Prop.('width_BMP'),zi*a); % the 0.1 ratio is empirical. There needs to be a better way
    El.dim(eX,2) = feval(Prop.(bz),zi*a); % depth
    
    eX = eBMza; bz = 'thick_BMp';    
    zi = mean([Nd.Z(El.Nd1(eX)),Nd.Z(El.Nd2(eX))], 2);
    El.dim(eX,1) = 0.333*feval(Prop.('width_BMA'),zi*a); % the 0.1 ratio is empirical. There needs to be a better way
    El.dim(eX,2) = 0.333*feval(Prop.(bz),zi*a); % depth
  
    
%     mthick_BMz = thickness_BMz(Nd, El, eX, MP);
    
    eX = eOHC;  bz = 'diam_OHC';                          
    zi = mean([Nd.Z(El.Nd1(eX)),Nd.Z(El.Nd2(eX))], 2);     
    El.dim(eX,1) = feval(Prop.(bz),zi*a);            
    El.dim(eX,2) = feval(Prop.(bz),zi*a);                   

    eX = eIPC; bz = 'thick_IPC';
    zi = mean([Nd.Z(El.Nd1(eX)),Nd.Z(El.Nd2(eX))], 2);    
    El.dim(eX,1) = feval(Prop.(bz),zi*a);             
    El.dim(eX,2) = feval(Prop.(bz),zi*a);   

    eX = eOPC; bz = 'thick_OPC';
    zi = mean([Nd.Z(El.Nd1(eX)),Nd.Z(El.Nd2(eX))], 2);    
    El.dim(eX,1) = feval(Prop.(bz),zi*a);        
    El.dim(eX,2) = feval(Prop.(bz),zi*a);   
    
    eX = eDCb;  bz = 'diam_DCb';                          
    zi = mean([Nd.Z(El.Nd1(eX)),Nd.Z(El.Nd2(eX))], 2);     
    El.dim(eX,1) = feval(Prop.(bz),zi*a);              
    El.dim(eX,2) = feval(Prop.(bz),zi*a);         
    
    eX = eDCp;  bz = 'diam_DCp';                          
    zi = mean([Nd.Z(El.Nd1(eX)),Nd.Z(El.Nd2(eX))], 2);    
    El.dim(eX,1) = feval(Prop.(bz),zi*a);                
    El.dim(eX,2) = feval(Prop.(bz),zi*a);          
   
    eX = eDCz; bz = 'diam_DCz';
    zi = mean([Nd.Z(El.Nd1(eX)),Nd.Z(El.Nd2(eX))], 2);
    El.dim(eX,1) = feval(Prop.('height_TOC'),zi*a).*(1 - feval(Prop.('rOHC_DC'),zi*a));
    El.dim(eX,2) = 3*feval(Prop.('diam_DCb'),zi*a);
    
    eX = eDCr;
    zi = mean([Nd.Z(El.Nd1(eX)),Nd.Z(El.Nd2(eX))], 2);
    El.dim(eX,1) = feval(Prop.(bz),zi*a);
    El.dim(eX,2) = feval(Prop.(bz),zi*a);
    
    eX = eRLx1; bz = 'thick_RL';
    zi = Nd.Z(El.Nd2(eX));                              % Nd1 to Nd2, MPN    
    El.dim(eX,1) = MP.dZ;
    El.dim(eX,2) = feval(Prop.(bz),zi*a);   
    
    eX = eRLp; bz = 'thick_RL';
    zi = Nd.Z(El.Nd2(eX));                              % Nd1 to Nd2, MPN
    El.dim(eX,1) = MP.dZ;
    El.dim(eX,2) = feval(Prop.(bz),zi*a);   
    
    eX = eRLpz; bz = 'thick_RL';
    zi = Nd.Z(El.Nd2(eX));                              % Nd1 to Nd2, MPN
    El.dim(eX,1) = feval(Prop.('diam_DCpr'),zi*a);
    El.dim(eX,2) = feval(Prop.(bz),zi*a);   

    eX = eRLz; bz = 'thick_RL';
    zi = mean([Nd.Z(El.Nd1(eX)),Nd.Z(El.Nd2(eX))], 2);    
    El.dim(eX,1) = 0.5*feval(Prop.('width_OHC_RL'),zi*a); % the 0.5 ratio is empirical. There needs to be a better way
    El.dim(eX,2) = feval(Prop.(bz),zi*a);
    
    eX = eTMx1; bz = 'thick_TMa';   % TM attachement thickness
    zi = mean([Nd.Z(El.Nd1(eX)),Nd.Z(El.Nd2(eX))], 2);
    El.dim(eX,1) = MP.dZ;    
    El.dim(eX,2) = feval(Prop.(bz),zi*a);   
    % attachement region has half the thickness of the body
    
    eX = eTMx2; bz = 'thick_TMb';   % TM body thickness
    zi = mean([Nd.Z(El.Nd1(eX)),Nd.Z(El.Nd2(eX))], 2);
    El.dim(eX,1) = MP.dZ;
    El.dim(eX,2) = feval(Prop.(bz),zi*a);   

    eX = eTMz0; bz = 'thick_TMb';                       % Longitudianl TM attachment
    zi = mean([Nd.Z(El.Nd1(eX)),Nd.Z(El.Nd2(eX))], 2);
    El.dim(eX,1) = 0.33*feval(Prop.('width_TM'),zi*a); % the 0.33 ratio is empirical. There needs to be a better way
    El.dim(eX,2) = feval(Prop.(bz),zi*a);   

    eX = eTMz1; bz = 'thick_TMb';                       % Longitudianl TM at body
    zi = mean([Nd.Z(El.Nd1(eX)),Nd.Z(El.Nd2(eX))], 2);
    El.dim(eX,1) = 0.66*feval(Prop.('width_TM'),zi*a); % the 0.66 ratio is empirical. There needs to be a better way
    El.dim(eX,2) = feval(Prop.(bz),zi*a);   
    
    
%N%%N%%N%%N%%N%%N%%N%%N%%N%%N%%N%%N%%N%%N%%N%%N%%N%%N%%N%%N%%N%%N%%N%%N%%N%%N%%N%%N%%N%%N%%N%%N%%N%
%N%%N%%N%%N%%N%%N%%N%%N%%N%%N%%N%%N%%N%%N%%N%%N%%N%%N%%N%%N%%N%%N%%N%%N%%N%%N%%N%%N%%N%%N%%N%%N%%N%
%N%%N%%N%%N%%N%%N%%N%%N%%N%%N%%N%%N%%N%%N%%N%%N%%N%%N%%N%%N%%N%%N%%N%%N%%N%%N%%N%%N%%N%%N%%N%%N%%N%
    eX = eOHB;    
    zi = mean([Nd.Z(El.Nd1(eX)),Nd.Z(El.Nd2(eX))], 2);
    bz = 'k_OHB';   kHB = feval(Prop.(bz),zi*a);   
    bz = 'height_HB';   hHB = feval(Prop.(bz),zi*a);   
    YM_HB = 10e6;
    a1 = 0.9;
        
    wHB = 10.0;  % HB width [um]
    El.dim(eX,1) = wHB;
    El.dim(eX,2) = hHB.*( (4*a1*kHB/YM_HB/wHB).^(1/3));        % Cantilever beam stiffness kHB = 3*YM*Iz/H^3, Iz = (w*t^3)/12 --> t = H*(4*kHB/YM/w)^1/3
    
    eX = exHB;
    El.dim(eX,1) = wHB;
    El.dim(eX,2) = 0.5;
    
    eX = eANK;
    El.dim(eX,1) = wHB;
    El.dim(eX,2) = (1-a1)*kHB.*El.oL(eX)/YM_HB/wHB;            % Bar stiffness kANK = YM*w*t/L, where L = El.oL --> t = kANK*L/YM/w
    
%N%%N%%N%%N%%N%%N%%N%%N%%N%%N%%N%%N%%N%%N%%N%%N%%N%%N%%N%%N%%N%%N%%N%%N%%N%%N%%N%%N%%N%%N%%N%%N%%N%
%N%%N%%N%%N%%N%%N%%N%%N%%N%%N%%N%%N%%N%%N%%N%%N%%N%%N%%N%%N%%N%%N%%N%%N%%N%%N%%N%%N%%N%%N%%N%%N%%N%
%N%%N%%N%%N%%N%%N%%N%%N%%N%%N%%N%%N%%N%%N%%N%%N%%N%%N%%N%%N%%N%%N%%N%%N%%N%%N%%N%%N%%N%%N%%N%%N%%N%    
    
    eX = eTCz; bz = 'thick_TCz';
    zi = mean([Nd.Z(El.Nd1(eX)),Nd.Z(El.Nd2(eX))], 2);
    El.dim(eX,1) = 0.5*feval(Prop.('width_OHC_RL'),zi*a); % the 0.5 ratio is empirical. There needs to be a better way
    El.dim(eX,2) = feval(Prop.(bz),zi*a);
    
end % of function set_dimension_gradient()
