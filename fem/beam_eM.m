function eM = beam_eM(El,ie,type)

    if ~exist('type','var')
        type = 'lumped';
    end
    % Note that these are inconsistent with K matrix
    % but, it would affect little to the dynamics
    A = El.A(ie);
    L = El.L(ie);
    a = 0.5*El.L(ie);
    rho = El.rho(ie);
    rx2 = 2*El.Iz(ie)/A;
    qq=rho*A*a/105;
    % handle transverse matrix first     

    if strcmp(type,'consistent')
        eM = qq*[   70, 0, 0, 0, 0, 0,          35, 0, 0, 0, 0, 0;
                    0, 78, 0, 0, 0, 22*a,       0, 27, 0, 0, 0, -13*a;
                    0, 0, 78, 0, -22*a, 0,      0, 0, 27, 0, 13*a, 0;
                    0, 0, 0, 70*rx2, 0, 0,      0, 0, 0, 35*rx2, 0, 0;
                    0, 0, -22*a, 0, 8*a^2, 0,   0, 0, -13*a, 0, -6*a^2, 0;
                    0, 22*a, 0, 0, 0, 8*a^2,    0, 13*a, 0, 0, 0, -6*a^2;
    
                    35, 0, 0, 0, 0, 0,          70, 0, 0, 0, 0, 0;
                    0, 27, 0, 0, 0, 13*a,       0, 78, 0, 0, 0, -22*a;
                    0, 0, 27, 0, -13*a, 0,      0, 0, 78, 0, 22*a, 0;
                    0, 0, 0, 35*rx2, 0, 0,      0, 0, 0, 70*rx2, 0, 0;
                    0, 0, 13*a, 0, -6*a^2, 0,   0, 0, 22*a, 0, 8*a^2, 0;
                    0, -13*a, 0, 0, 0, -6*a^2,  0, -22*a, 0, 0, 0, 8*a^2;];
    elseif strcmp(type,'lumped')
        % lumped mass with zero rotary inertia
        eM = rho*A*L*diag([0.5,0.5,0.5,0,0,0,0.5,0.5,0.5,0,0,0]);
    else
        error('Mass matrix type not defined');
    end
    % lumped mass with zero rotary inertia
    if contains(El.name(ie),'TM')
        eM = rho*A*L*diag([0.5,0.5,0.5,0,0,0,0.5,0.5,0.5,0,0,0]);
    end

    R = rotation_matrix(El,ie);        
            
    eM = R'*eM*R;
    
end


function  R = rotation_matrix(El,ie)
    
    if abs(El.dir(ie,3)-1)<0.1,
        ez = [1, 0, 0];
        % note that most elements alligned either radial(x) or longitudianl(z) direction
    else
        ez = [0, 0, 1];
    end

    xv = El.dir(ie,:);

    yv = mycross(ez,xv);
    yv = yv/norm(yv);
    zv = mycross(xv,yv);
    zv = zv/norm(zv);

    T = [xv;yv;zv];
% %     R = zeros(9);
% %     for ii=1:4,
% %         R((ii-1)*3+1:(ii-1)*3+3,(ii-1)*3+1:(ii-1)*3+3) = T;
% %     end
    O33 = zeros(3);    
    R = [T, O33, O33, O33;  O33, T, O33, O33;  O33, O33, T, O33;   O33, O33, O33, T;];
end