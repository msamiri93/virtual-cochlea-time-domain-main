function El = get_vector(El, Nd, opt)
%   function [len, dir] = get_vector(El, Nd)
% -----------------------------------------------------------------------
% 
%   Program unit: get_vector
%   Description: get distance and direction vector between two nodes
%
%   First written: 2009.03.25
%
%   by Jong-Hoon Nam
% 

    for ie = 1:El.N,
        
        if El.type(ie)<2.5, % if it is a 1D element

            P1 = [ Nd.X(El.Nd1(ie)); Nd.Y(El.Nd1(ie)); Nd.Z(El.Nd1(ie)) ]; % coord of node 1
            P2 = [ Nd.X(El.Nd2(ie)); Nd.Y(El.Nd2(ie)); Nd.Z(El.Nd2(ie)) ]; % coord of node 2

            D1 = [ Nd.dx(El.Nd1(ie)); Nd.dy(El.Nd1(ie)); Nd.dz(El.Nd1(ie)) ]; % coord of node 1
            D2 = [ Nd.dx(El.Nd2(ie)); Nd.dy(El.Nd2(ie)); Nd.dz(El.Nd2(ie)) ]; % coord of node 2

            P1 = P1 + D1;
            P2 = P2 + D2;

            distance = norm(P2-P1);
            direction = (P2-P1)/norm(P2-P1);

            if opt == 0,
                El.oL(ie,1) = distance;
                El.odir(ie,:) = direction;
            end

            El.L(ie,1) = distance;
            El.dir(ie,:) = direction;
            
        end
    end
end