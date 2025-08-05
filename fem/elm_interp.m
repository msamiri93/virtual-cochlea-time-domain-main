function [psi, gpsi, hs] = elm_interp(type, xx, xi, eta)
% Evaluation the surface metric coefficient, $h_s$,
% and computation of the basis functions and their gradients,
% over a quadrilateral.

    if type == 4

        %%%% compute the basis functions
        psi = [((eta - 1)*(xi - 1))/4, -((eta - 1)*(xi + 1))/4, ...
            ((eta + 1)*(xi + 1))/4, -((eta + 1)*(xi - 1))/4];
        %%%% compute the xi derivatives of the basis functions
        dps = [eta/4 - 1/4, 1/4 - eta/4, eta/4 + 1/4, - eta/4 - 1/4];

        %%%% compute the eta derivatives of the basis functions
        pps = [xi/4 - 1/4, - xi/4 - 1/4, xi/4 + 1/4, 1/4 - xi/4];

    elseif type == 8
        
        %%%% compute the basis functions
        psi = [ -((eta - 1)*(xi - 1)*(eta + xi + 1))/4, ((eta - 1)*(xi + ...
            1)*(eta - xi + 1))/4, ((eta + 1)*(xi + 1)*(eta + xi - 1))/4, ...
            ((eta + 1)*(xi - 1)*(xi - eta + 1))/4, (xi^2/2 - 1/2)*(eta - 1)...
            , -(eta^2 - 1)*(xi/2 + 1/2), -(xi^2/2 - 1/2)*(eta + 1), ...
            (eta^2 - 1)*(xi/2 - 1/2)];
        %%%% compute the xi derivatives of the basis functions
        dps = [ -((eta + 2*xi)*(eta - 1))/4, ((eta - 2*xi)*(eta - 1))/4, ...
            ((eta + 2*xi)*(eta + 1))/4, -((eta - 2*xi)*(eta + 1))/4, ...
            xi*(eta - 1), 1/2 - eta^2/2, -xi*(eta + 1), eta^2/2 - 1/2];

        %%%% compute the eta derivatives of the basis functions
        pps = [ -((2*eta + xi)*(xi - 1))/4, ((xi + 1)*(2*eta - xi))/4, ...
            ((2*eta + xi)*(xi + 1))/4, -((xi - 1)*(2*eta - xi))/4, xi^2/2 ...
            - 1/2, -eta*(xi + 1), 1/2 - xi^2/2, eta*(xi - 1)];
        
    elseif type == 9
        
        %%%% compute the basis functions
        psi = [ (eta*xi*(eta - 1)*(xi - 1))/4, (eta*xi*(eta - 1)*(xi + 1))/4 ...
            , (eta*xi*(eta + 1)*(xi + 1))/4, (eta*xi*(eta + 1)*(xi - 1))/4, ...
            -(eta*(xi^2 - 1)*(eta - 1))/2, -(xi*(eta^2 - 1)*(xi + 1))/2, ...
            -(eta*(xi^2 - 1)*(eta + 1))/2, -(xi*(eta^2 - 1)*(xi - 1))/2, ...
            (eta^2 - 1)*(xi^2 - 1)];
        
        %%%% compute the xi derivatives of the basis functions
        dps = [ (eta*(2*xi - 1)*(eta - 1))/4, (eta*(2*xi + 1)*(eta - 1))/4, ...
            (eta*(2*xi + 1)*(eta + 1))/4, (eta*(2*xi - 1)*(eta + 1))/4, ...
            -eta*xi*(eta - 1), -((eta^2 - 1)*(2*xi + 1))/2, -eta*xi*(eta + 1), ...
            -((eta^2 - 1)*(2*xi - 1))/2, 2*xi*(eta^2 - 1)];
        
        %%%% compute the eta derivatives of the basis functions
        pps = [ (xi*(2*eta - 1)*(xi - 1))/4, (xi*(2*eta - 1)*(xi + 1))/4, ...
            (xi*(2*eta + 1)*(xi + 1))/4, (xi*(2*eta + 1)*(xi - 1))/4, ...
            -((2*eta - 1)*(xi^2 - 1))/2, -eta*xi*(xi + 1), ...
            -((2*eta + 1)*(xi^2 - 1))/2, -eta*xi*(xi - 1), 2*eta*(xi^2 - 1)];
                
    end

    %%%% compute the xi and eta derivatives of x
    DxDxi = dps*xx;
    DxDet = pps*xx;      

    %%%% The system is solved by Cramer's rule    
    A = [DxDxi; DxDet];    % a 2-by-2 matrix
    bb = [dps; pps];       % a 2-by-8 matrix

    %%%% compute the surface metric hs        
    hs = abs(det(A)); 

    gpsi = (inv(A)*bb)';        % a 8-by-2 vector
    
end