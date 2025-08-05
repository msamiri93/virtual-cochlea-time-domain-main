function [Ke, R] = beam_eK(El, ie)
   
    R = rotation_matrix(El,ie);

    if El.type(ie) == 0,
        Ke = timochenko_beam(El, ie, R);
    elseif El.type(ie) == 0.1
        Ke = euler_beam(El, ie, R);
    else        
%         Ke = euler_beam_hybrid(El, ie, R);
    end
end
    
function  R = rotation_matrix(El,ie)
    
    if abs(El.dir(ie,3)-1)<0.1,
        ez = [1, 0, 0];
        % note that most elements alligned either radial(x) direction
        % or longitudianl(z) direction
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
    

function Ke = timochenko_beam(El, ie, R)
        
        SM = El.SM(ie);
        YM = El.YM(ie);
        A = El.A(ie);
        Iz = El.Iz(ie);
        Iy = El.Iy(ie);
        Jx = Iz + Iy;

        L = El.oL(ie);

        A1=SM*A/L; 
        B1=-SM*A/2; 
        C1=SM*A*L/4; 
        D1=YM*Iy/L; 
        D2=YM*Iz/L;
        E1=YM*A/L; 
        F1 = SM*Jx/L;

        Ke = [E1, 0, 0, 0, 0, 0,        -E1, 0,  0, 0, 0, 0;
              0,  A1, 0, 0, 0, B1,        0,-A1, 0, 0, 0, B1;
              0, 0, A1, 0, B1, 0,         0, 0,-A1, 0, B1, 0;
              0, 0, 0, F1, 0, 0,          0, 0,  0,-F1, 0, 0;
              0, 0, B1, 0, C1+D1, 0,      0, 0,-B1, 0, C1-D1, 0;
              0, B1, 0, 0, 0, C1+D2,      0,-B1, 0, 0, 0, C1-D2;


              -E1, 0, 0, 0, 0, 0,        E1, 0, 0, 0, 0, 0;
              0, -A1, 0, 0, 0, -B1,       0,A1, 0, 0, 0 -B1;
              0, 0, -A1, 0, -B1, 0,       0, 0,A1, 0,-B1, 0;
              0, 0, 0, -F1, 0, 0,         0, 0, 0, F1, 0, 0;
              0, 0, B1, 0, C1-D1, 0,      0, 0,-B1, 0, C1+D1, 0;
              0, B1, 0, 0, 0, C1-D2,      0,-B1, 0, 0, 0, C1+D2 ]; 
          
%         if isfield(El,'hinge'),
%             if sum(El.hinge(ie,:))>0,
%                 Ke = install_hinge(Ke, El, ie);
%             end
%         end

        Ke = R'*Ke*R;
end
    
function Ke = euler_beam(El, ie, R) 
      
    YM = El.YM(ie);
    SM = El.SM(ie);
    A = El.A(ie);
    Iz = El.Iz(ie);
    Iy = El.Iy(ie);
    Jx = Iz + Iy;
   
    L = El.oL(ie);

    A1 = A*YM/L;
    B1 = 12*YM*Iz/L^3;
    C1 = 12*YM*Iy/L^3;
    D1 = SM*Jx/L;
    E1 = 4*YM*Iy/L;
    F1 = 4*YM*Iz/L;
    G1 = 6*YM*Iy/L^2;
    H1 = 6*YM*Iz/L^2;


    Ke = [ A1,   0,   0,   0,   0,   0,    -A1,   0,   0,   0,   0,   0;
            0,  B1,   0,   0,   0,  H1,      0, -B1,   0,   0,   0,  H1;
            0,   0,  C1,   0, -G1,   0,      0,   0, -C1,   0, -G1,   0;
            0,   0,   0,  D1,   0,   0,      0,   0,   0, -D1,   0,   0;
            0,   0, -G1,   0,  E1,   0,      0,   0,  G1,   0, E1/2,  0;
            0,  H1,   0,   0,   0,  F1,      0, -H1,   0,   0,   0, F1/2;


          -A1,   0,   0,   0,   0,   0,     A1,   0,   0,   0,   0,   0;
            0, -B1,   0,   0,   0, -H1,      0,  B1,   0,   0,   0, -H1;
            0,   0, -C1,   0,  G1,   0,      0,   0,  C1,   0,  G1,   0;
            0,   0,   0, -D1,   0,   0,      0,   0,   0,  D1,   0,   0;
            0,   0, -G1,   0, E1/2,  0,      0,   0,  G1,   0,  E1,   0;
            0,  H1,   0,   0,   0, F1/2,     0, -H1,   0,   0,   0,  F1 ]; 
    
%     if isfield(El,'hinge'),
%             if sum(El.hinge(ie,:))>0,
%                 Ke = install_hinge(Ke, El, ie);
%             end
%     end

    Ke = R'*Ke*R;
    
end

function Ke = euler_beam_hybrid(El, ie, R)

    global coef;

    YM = El.YM(ie);
    SM = El.SM(ie);
    A = El.A(ie);
    Iz = El.Iz(ie);
    Iy = El.Iy(ie);
    Jx = Iz + Iy;
   
    L = El.oL(ie);

    A1 = A*YM/L;
    B1 = 12*YM*Iz/L^3;
    C1 = 12*YM*Iy/L^3;
    D1 = SM*Jx/L;
    E1 = 4*YM*Iy/L;
    F1 = 4*YM*Iz/L;
    G1 = 6*YM*Iy/L^2;
    H1 = 6*YM*Iz/L^2;


    Ke = [ A1,   0,   0,   0,   0,   0,    -A1,   0,   0,   0,   0,   0;
            0,  B1,   0,   0,   0,  H1,      0, -B1,   0,   0,   0,  H1;
            0,   0,  C1,   0, -G1,   0,      0,   0, -C1,   0, -G1,   0;
            0,   0,   0,  D1,   0,   0,      0,   0,   0, -D1,   0,   0;
            0,   0, -G1,   0,  E1,   0,      0,   0,  G1,   0, E1/2,  0;
            0,  H1,   0,   0,   0,  F1,      0, -H1,   0,   0,   0, F1/2;


          -A1,   0,   0,   0,   0,   0,     A1,   0,   0,   0,   0,   0;
            0, -B1,   0,   0,   0, -H1,      0,  B1,   0,   0,   0, -H1;
            0,   0, -C1,   0,  G1,   0,      0,   0,  C1,   0,  G1,   0;
            0,   0,   0, -D1,   0,   0,      0,   0,   0,  D1,   0,   0;
            0,   0, -G1,   0, E1/2,  0,      0,   0,  G1,   0,  E1,   0;
            0,  H1,   0,   0,   0, F1/2,     0, -H1,   0,   0,   0,  F1 ]; 
    
    if El.type > 0, n = 2; else, n = 1; end
    El.hinge(ie,n) = 1;
    Ke_h = install_hinge(Ke, El, ie);

    a = coef.hinge_ratio;
    Ke = (1-a)*Ke + a*Ke_h;

    Ke = R'*Ke*R;

end
    
    
 function Ke = install_hinge(Ke, El, ie)

     if El.hinge(ie,1)==1,
        jdx = ones(1,12); jdx(4:6)=0; a1 = jdx==1; a2 = jdx==0;
        K11 = Ke(a1,a1);
        K12 = Ke(a1,a2);
        K21 = Ke(a2,a1);
        K22 = Ke(a2,a2);
        Keq = K11 - K12*(K22\K21);
        
        Ke = zeros(12,12);
        Ke(a1,a1) = Keq;        
     elseif El.hinge(ie,2)==1,
        jdx = ones(1,12); jdx(10:12)=0; a1 = jdx==1; a2 = jdx==0;
        K11 = Ke(a1,a1);
        K12 = Ke(a1,a2);
        K21 = Ke(a2,a1);
        K22 = Ke(a2,a2);
        Keq = K11 - K12*(K22\K21);
        
        Ke = zeros(12,12);
        Ke(a1,a1) = Keq;         
     else
         fprintf(1,'\n Element''s hinged condition is inappropriate. No hinge installed.\n');
     end        
 end