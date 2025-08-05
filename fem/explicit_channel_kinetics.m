function po = explicit_channel_kinetics(time, MP, HB, po, xHB, xa)
        
        options = odeset('Jacobian',jFunc(xHB, xa, HB),'RelTol',1e-3,'AbsTol',1e-6,'Stats','off');
        [t,y] = ode23(@(t,y) channel_2state_t(t, y, xHB, HB, xa),[time time+MP.dt],po,options);

        po = y(end,:);
end

function dpodt = channel_2state_t(t, po, xx, HB, xa)

    po = reshape(po,1,[]);
    xx = reshape(xx,1,[]);
    xa = reshape(xa,1,[]);

    dE = [HB.z].*(xx - xa - [HB.Xo]);
    kBT = 4e-3;
    kCO = [HB.kF].*exp(0.5*dE/kBT);      % rate from closed to open state
    kOC = [HB.kR].*exp(-0.5*dE/kBT);     % rate from open to closed state
    dpodt = reshape(kCO.*(1-po)-kOC.*po,[],1);

end

function dfdpo = jFunc(xx, xa, HB)

    nn = numel(xx);
    xx = reshape(xx,1,[]);
    xa = reshape(xa,1,[]);  
    dE = [HB.z].*(xx - xa - [HB.Xo]);
    kBT = 4e-3;
    kCO = [HB.kF].*exp(0.5*dE/kBT);      % rate from closed to open state
    kOC = [HB.kR].*exp(-0.5*dE/kBT);     % rate from open to closed state
    dfdpo = spdiags(-(kCO + kOC),0,nn,nn);

end