function [xa, po, dt_s] = explicit_channel_kinetics_r02(time, dt_s, MP, HB, po, xHB, xa)
        
        options = odeset('Jacobian',@(t,y) jFunc(t, y, xHB, HB),'InitialStep',dt_s,'RelTol',1e-4,'AbsTol',1e-8,'Stats','off');
        y0 = [xa;po];
        [t,y] = ode23(@(t,y) channel_2state_t(t, y, xHB, HB),[time time+MP.dt],y0,options);
        dt_s = t(end) - t(end-1);
        y = transpose(reshape(y(end,:),[],2));
        xa = y(1,:);
        po = y(2,:);
end

function dydt = channel_2state_t(t, y, xx, HB)

    y = reshape(y,[],2);
    xa = reshape(y(:,1),1,[]);
    po = reshape(y(:,2),1,[]);
    xx = reshape(xx,1,[]);
    

    dE = [HB.z].*(xx - xa - [HB.Xo]);
    kBT = 4e-3;
    kCO = [HB.kF].*exp(0.5*dE/kBT);      % rate from closed to open state
    kOC = [HB.kR].*exp(-0.5*dE/kBT);     % rate from open to closed state
    dxadt = reshape([HB.kA].*([HB.kG].*[HB.gamma].*(xx - xa) -[HB.kES].*xa),[],1);
    dpodt = reshape(kCO.*(1-po)-kOC.*po,[],1);
    dydt = [dxadt;dpodt];
end

function JJ = jFunc(t,y, xx, HB)

    nn = numel(xx);
    y = reshape(y,[],2);
    xa = reshape(y(:,1),1,[]);
    po = reshape(y(:,2),1,[]);

    xx = reshape(xx,1,[]);

    dE = [HB.z].*(xx - xa - [HB.Xo]);
    kBT = 4e-3;
    kCO = [HB.kF].*exp(0.5*dE/kBT);      % rate from closed to open state
    kOC = [HB.kR].*exp(-0.5*dE/kBT);     % rate from open to closed state

    % Jacobian elements
    J11 = -[HB.kA] .* ( [HB.kG].*[HB.gamma] + [HB.kES]);  % d(dxadt)/dxa
    J12 = zeros(1,nn);                  % d(dxadt)/dpo
    J21 = -0.5 * [HB.z] ./ kBT .* (kCO .* (1 - po) + kOC .* po); % d(dpodt)/dxa
    J22 = -kCO - kOC;         % d(dpodt)/dpo

    JJ = [spdiags(J11,0,nn,nn),spdiags(J12,0,nn,nn);
          spdiags(J21,0,nn,nn),spdiags(J22,0,nn,nn)];

end