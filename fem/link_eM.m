function eM = link_eM(El, ie,type)

    if ~exist('type','var')
        type = 'lumped';
    end

    rA = El.rho(ie)*El.A(ie);        
    hL = 0.5*El.L(ie);    
    if strcmp(type,'consistent')
        eM = [  70, 0, 0, 35, 0, 0;
                0, 78, 0, 0, 27, 0;
                0, 0, 78, 0, 0, 27;
                35, 0, 0, 70, 0, 0;
                0, 27, 0, 0, 78, 0;
                0, 0, 27, 0, 0, 78;];
        eM = rA*hL/105*eM;

    elseif strcmp(type,'lumped')
        eM = rA*hL*diag([1,1,1,1,1,1]);
    else
        error('Mass matrix type not defined');
    end   

    R = rotation_matrix(El, ie);            
    eM = R'*eM*R;   

end


function R = rotation_matrix(El, ie)

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

    T33 = [xv;yv;zv];
    O33 = zeros(3);
    
    R = [T33, O33; O33, T33];

end