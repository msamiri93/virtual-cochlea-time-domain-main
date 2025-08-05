function plot_cross_section(ax,MP,Nd,El,U,options)
% loc: location
% nf: # reference of the frequency

    arguments
        ax
        MP
        Nd
        El
        U
        options.loc = mean(MP.xx);
        options.freq = 0;
        options.cycle = -1;
        options.pfact = 10;
    end

    loc = options.loc;
    freq = options.freq;

    NDOF = 6;

    pfact = options.pfact;
    % if freq == 0
    %     pfact = 0;
    % else
    %     pfact = options.pfact;
    % end

    nt = size(U,2);
    nz = knnsearch(MP.xx(:),loc);
    nn = ((1:Nd.Nc)-1)*Nd.Nr + nz;

    xx = Nd.X(nn);
    yy = Nd.Y(nn);
    zz = Nd.Z(nn);

    if isfield(MP,'dof')
        udof = MP.dof(2).udof;
        uu = reshape(U(udof,:),NDOF,Nd.N,nt);
    else
        uu = reshape(U(:,nf),NDOF,Nd.N);
    end
    dx = uu(1,nn,:);
    dy = uu(2,nn,:);
    dz = uu(3,nn,:);
    
    ps = pfact/max([abs(dx(:));abs(dy(:))]);

    hold(ax,'on');

    n_sec = find(abs(Nd.Z + MP.loc*1e3 - loc*1e3) < 15);
    e_idx = find(ismember(El.Nd1,nn));

    opt_anim = 0;
    if options.cycle > 0
        np = 41;
        phi = linspace(-pi,pi,np);% + pi/np*rand(1,np);  
        if options.cycle > 1
            opt_anim = 1;
        end
    elseif options.cycle == -1
        nAF = find(strcmp(Nd.name,'AF'));
        phi = pi - angle(U(nAF(nz),nf));
    else
        phi = options.cycle;        
    end

    icount = 0;

    if opt_anim == 1
        file_name = './houtput/anim.gif';
        if exist(file_name,'file'), delete(file_name); end        
    end
    for si = phi
        ti = round(interp1([-pi pi],[1 nt],si));
        icount = icount + 1;
        for ie=1:El.N
            if contains(El.name(ie),'z')
                continue
            end
            if ismember(ie,e_idx)
            nds = [El.Nd1(ie), El.Nd2(ie)];               
                %undeformed shadow
                [R,T] = rotation_matrix(El,ie);
                % d = uu(:,nds)*exp(1i*si);
                d = squeeze(uu(:,nds,ti));
                nl = 20;
                ang = linspace(angle(d(2,1)),angle(d(2,2)),nl);
                coord = [Nd.X(nds);Nd.Y(nds);Nd.Z(nds)];
                d = R*d(:); % move displacement values to local frame
                L = El.L(ie);
                x = linspace(0,El.L(ie),nl);
                xx_l = zeros(3,nl);
                dy_l = zeros(1,nl);
                dx_l = zeros(1,nl);
                dz_l = zeros(1,nl);
                O = zeros(1,nl);
                for ii = 1:numel(x)
                    % this is in x-y plane
                    N = [(L-x(ii))/L, x(ii)/L]; %[u1,u2]
                    xx_l(:,ii) = coord*N.';
                    dx_l(ii) = N*d([1;7]);
                    N = [1-3*x(ii)^2/L^2+2*x(ii)^3/L^3, x(ii)-2*x(ii)^2/L+x(ii)^3/L^2,...
                         3*x(ii)^2/L^2-2*x(ii)^3/L^3,-x(ii)^2/L+x(ii)^3/L^2]; %[v1,t1,v2,t2] local
                    dy_l(ii) = real(N*d([2;6;8;12]));
                    dz_l(ii) = real(N*d([3;5;9;11]));
                end
                dd = T.'*ps*[dx_l;dy_l;dz_l];
                % plot(ax,xx_l(1,:) + dd(1,:), xx_l(2,:) + dd(2,:),'Color',el_color,'LineWidth', 1.5);
                
                el_color = define_element_color(ie,El,ang,phi,options.cycle);

                if options.cycle == 1
                    plot3(ax,xx_l(3,:) + real(dd(3,:)), xx_l(1,:) + real(dd(1,:)) ...
                        , xx_l(2,:) + real(dd(2,:)),'Color',el_color,'LineWidth', 1.5);
                else
                    lw = 2;
                    if contains(El.name(ie),'OHC')
                        lw = 3*lw;
                    end
                    plot(ax,xx_l(1,:), xx_l(2,:),'Color',[0.6,0.6,0.6,0.4],'LineWidth',2);
                    patch(ax,[xx_l(1,:) + real(dd(1,:)),NaN], [xx_l(2,:) + real(dd(2,:)),NaN]...
                        , 1,'FaceVertexCData',[el_color;[NaN,NaN,NaN]],...
                        'EdgeColor','interp','LineWidth',lw);
                end

            end
        end

        if options.cycle == 1
            % plot3(ax,zz + ps*real(dz*exp(1i*si)), xx + ps*real(dx*exp(1i*si)), ...
            %     yy + ps*real(dy*exp(1i*si)), 'o','Color','none','MarkerFaceCol', ...
            %     [0, 0.4470, 0.7410],'MarkerSize',4);
            plot3(ax,zz + ps*real(squeeze(dz(:,:,ti))), xx + ps*real(squeeze(dx(:,:,ti))), ...
            yy + ps*real(squeeze(dy(:,:,ti))), 'o','Color','none','MarkerFaceCol', ...
            [0, 0.4470, 0.7410],'MarkerSize',4);
        end

        daspect(ax,[1,1,1])
        
        if opt_anim == 1
            xlim(ax,[-200 350])
            ylim(ax,[-50 150])
            xticks(ax,[])
            yticks(ax,[])
            exportgraphics(ax,file_name,Append=true);
            cla(ax);
        end
    end
    
end

function el_color = define_element_color(ie,El,ang,phi,cycle)

    cmap = hsv(256);
    if cycle == 1
        if contains(El.name(ie),'BM') || contains(El.name(ie),'TM') 
            el_color = [0.0,0.447,0.741,0.2];
        elseif contains(El.name(ie),'OHC')   
            el_color = [0.929,0.694,0.125,0.2];
        elseif contains(El.name(ie),'DC')   
            el_color = [0.850,0.325,0.098,0.2];
        elseif contains(El.name(ie),'OHB')   
            el_color = [0.494,0.184,0.556,0.2];
        else
            el_color = [0.6,0.6,0.6,0.4];
        end
    else
        el_color = cmap(round(interp1([-pi,pi],[1,256],abs(ang),'linear')),:);
    end

end

function  [R,T] = rotation_matrix(El,ie)
    
    if abs(El.dir(ie,3)-1)<0.1
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