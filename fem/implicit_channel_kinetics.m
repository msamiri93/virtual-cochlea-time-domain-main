function [xa, po, dxadt, dpodt, dt_s] = implicit_channel_kinetics(time, dt_s, MP, HB, po, xHB, xa, dxadt, dpodt)
        
        options = odeset('Jacobian',@(t, y, yp) jFunc(t, y, yp, xHB, HB),'InitialStep',dt_s,'RelTol',1e-6,'AbsTol',1e-12,'Stats','off');
        y0 = [xa;po];
        y0p = [dxadt;dpodt];
        [t,y] = ode15i(@(t,y,yp) channel_2state_t_i(t, y, yp, xHB, HB),[time time+MP.dt],y0,y0p,options);
        if numel(t) > 1
            dt_s = t(end) - t(end-1);
        end
        yp = channel_2state_t([], y(end,:), xHB, HB);

        y = transpose(reshape(y(end,:),[],2));
        xa = y(1,:);
        po = y(2,:);

        yp = reshape(yp,[],2);
        dxadt = yp(:,1);
        dpodt = yp(:,2);   
end

function dy = channel_2state_t_i(~, y, yp, xx, HB)

    y = reshape(y,[],2);
    xa = reshape(y(:,1),1,[]);
    po = reshape(y(:,2),1,[]);          
    xx = reshape(xx,1,[]);

    yp = reshape(yp,[],2);
    dxadt = yp(:,1);
    dpodt = yp(:,2);

    dE = [HB.z].*(xx - xa - [HB.Xo]);
    kBT = 4e-3;
    kCO = [HB.kF].*exp(0.5*dE/kBT);      % rate from closed to open state
    kOC = [HB.kR].*exp(-0.5*dE/kBT);     % rate from open to closed state
    dxa = reshape([HB.kA].*([HB.kG].*[HB.gamma].*(xx - xa) -[HB.kES].*xa),[],1) - dxadt;
    dpo = reshape(kCO.*(1-po)-kOC.*po,[],1) - dpodt;
    dy = [dxa;dpo];
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

function [JJ,JP] = jFunc(~, y, ~, xx, HB)

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
    J13 = -ones(1,nn);                  % d(dxa)/dxadt
    J14 = zeros(1,nn);                  % d(dxa)/dpodt

    J21 = -0.5 * [HB.z] ./ kBT .* (kCO .* (1 - po) + kOC .* po); % d(dpodt)/dxa
    J22 = -kCO - kOC;         % d(dpodt)/dpo
    J23 = zeros(1,nn);        % d(dpo)/dxadt
    J24 = -ones(1,nn);         % d(dpo)/dpodt

    JJ = [spdiags(J11,0,nn,nn),spdiags(J12,0,nn,nn);
          spdiags(J21,0,nn,nn),spdiags(J22,0,nn,nn)];
    JP = [spdiags(J13,0,nn,nn),spdiags(J14,0,nn,nn);
          spdiags(J23,0,nn,nn),spdiags(J24,0,nn,nn)];

end