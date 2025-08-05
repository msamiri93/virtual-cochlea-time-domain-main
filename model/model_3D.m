function [Nd, El, MP] = model_3D(MP, dim, opt_lgradient, opt_MP_set)
% #
% # function [Nd, El, MP] = model_3D(MP, dim, opt_lgradient)
% #
% # Units are in pN, micrometer and Pa (N/m^2 or pN/um^2)
% #
% #
% # loc = location of partition 'b', 'a'
% # dim = 3, 2 (3dimension, 2 dimension)
% # opt_lgradient = optioin of longitudinal gradient 0 = no longitudinal
% #  gradient, 1 = gradient exists along longitudinal axis
% #
% # First written: 2009.06.xx
% # Last modified: 2021.07.xx
% #
% # by: Jong-Hoon Nam at University of Rochester
% #

MP.opt_lgradient = opt_lgradient;
MP.NDOF = 6;

if exist('dim','var')
    MP.dim = dim;   % 2D model or 3D model
else
    MP.dim = 3;
end

if isnumeric(MP.loc)
    if MP.loc <= 6
        location = 'b';
    else
        location = 'a';
    end
else
    fprintf(1,'\nError in model_3D, ''loc'' should be either numeric, ''a'' or ''b'' \n');
    MP.loc = 2;
    location = 'b';
end


if opt_MP_set == 0  
    MP = model_geometry(MP, location, opt_lgradient);    
    [Nd, MP] = model_nodes_og(MP);
    El = model_elements_og(Nd, MP);
    El.hinge = false(El.N,2);
    Nd.Master_Slave = zeros(sum(El.hinge(:)),2);
    % [Nd,El] = define_hinge_nodes(Nd,El);
    El = model_moduli5d_og(El, Nd, MP);    
else % This is a new set of geometry and material props updated on Apr 2021    
    [MP, Nd] = model_nodes_24(MP, opt_lgradient);
    El = model_elements_24(Nd, MP);
    El.hinge = false(El.N,2);
    Nd.Master_Slave = zeros(sum(El.hinge(:)),2);
    [Nd,El] = define_hinge_nodes(Nd,El);
    El = model_moduli5_24a(El, Nd, MP, opt_lgradient);   
end

Nd.Z = Nd.Z - 1e3*MP.loc;
if sum(El.YM == 123e6) > 0
    error('Undefined element YM in the model.');
end

if sum(El.dim < 0,'all')
    error('Negative element dimensions detected.');
end

% Nd = model_pArea(Nd, El);

MP.xOHC = 0;
MP.fOHC = 0;

end