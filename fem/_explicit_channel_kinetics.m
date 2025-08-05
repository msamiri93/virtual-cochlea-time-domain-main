function po = explicit_channel_kinetics(time, MP, HB, po, xHB, xa)

    % icount = 1;
    % t_i = time;
    % relTol = 1e-3;
    % % dt = MP.dt;
    % hmin = 16*eps(MP.dt);
    % dt = hmin;

    for si = 1:numel(po)

        icount = 1;
        t_i = time;
        relTol = 1e-4;
        % dt = MP.dt;
        hmin = 16*eps(MP.dt);
        dt = hmin;

        while true
            if icount > 1
                k1 = k4;
            else
                k1 = channel_2state_t(po(si), xHB(si), xa(si), HB(si));
            end
            k2 = channel_2state_t(po(si) + 0.5*dt*transpose(k1), xHB(si), xa(si), HB(si));
            k3 = channel_2state_t(po(si) + 0.75*dt*transpose(k2), xHB(si), xa(si), HB(si));
            z = po(si) + (2/9)*dt*transpose(k1) + (1/3)*dt*transpose(k2) + (4/9)*dt*transpose(k3);
            % z(z > 1) = 1; z(z < 0) = 0; % saturation correction
            k4 = channel_2state_t(z, xHB(si), xa(si), HB(si));
            po_est = po(si) + (7/24)*dt*transpose(k1) + (1/4)*dt*transpose(k2) + (1/3)*dt*transpose(k3) + (1/8)*dt*transpose(k4);
            % po_est(po_est > 1) = 1; po_est(po_est < 0) = 0; % saturation correction
            err = norm(z - po_est)/sqrt(mean(po.^2));
    
    
            
            if err > relTol            
                dt = max(hmin,dt* max(0.5,0.8*(relTol / err)^(1/3)));
                continue
            end
    
            if (t_i + dt) > (time + MP.dt)
                po(si) = transpose(interp1([t_i,t_i + dt],transpose([po(si),po_est]),time + MP.dt,'linear'));
                break;
            end    
    
            t_i = t_i + dt;
            icount = icount + 1;
            if err > 0
                % dt = dt*0.8*(relTol / err)^(1/3);
                dt = max(hmin,dt* min(2,0.8*(relTol / err)^(1/3)));
            else
                dt = 2*dt;
            end
            % dt = dt*0.8*(relTol / err)^(1/3);
            po(si) = po_est;
        end
    end
end

function dpodt = channel_2state_t(po, xx, xa, HB)

    po = reshape(po,1,[]);
    xx = reshape(xx,1,[]);
    xa = reshape(xa,1,[]);  
    dE = [HB.z].*(xx - xa - [HB.Xo]);
    kBT = 4e-3;
    kCO = [HB.kF].*exp(0.5*dE/kBT);      % rate from closed to open state
    kOC = [HB.kR].*exp(-0.5*dE/kBT);     % rate from open to closed state
    dpodt = kCO.*(1-po)-kOC.*po;

end