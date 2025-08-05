function eK = link_eK(El, ie)

    stiff = El.YM(ie)*El.A(ie)/El.oL(ie);
    eK = zeros(6);
    eK(1,1) = stiff;
    eK(4,4) = stiff;
    eK(1,4) = -1*stiff;
    eK(4,1) = -1*stiff;
    
    R = rotation_matrix(El, ie);
    
    eK = R'*eK*R;

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
