function Prop_fit = interpolate_properties(Prop_input,options)

    arguments
        Prop_input
        options.plot = true;
    end

    fnames = {Prop_input.name};

    if options.plot    
        f = figure("Name",'Geometric Properties');
        tabgp = uitabgroup(f,'Position',[.05 .05 .9 .9]);    
    end

    for i = 1:numel(fnames)        

        name = fnames{i};
        
        x = Prop_input(i).location;
        y = Prop_input(i).value;
        fit_type = Prop_input(i).fit_type;
        idx = ~isnan(y);
        x = x(idx);
        y = y(idx);

        if numel(unique(y)) == 1
            cfit = fit(x(:),y(:),'poly1');
            Prop_fit.(fnames{i}) = @(x) (feval(cfit,x));
            clear cfit;
            %warning('zero value parameter detected.');
        else
            % % % % % % if strcmp(name,'Y_OHC')
            % % % % % %     load('./geom/rOCC_OHC.mat','r','xi');
            % % % % % %     cfit = fit(xi(:),r(:),'smoothingspline');
            % % % % % %     Prop_fit.(fnames{i}) = @(x) ((mean(y)*feval(cfit,x)));
            % % % % % %     continue;
            % % % % % % end
    
            % exceptions check! this could be addressed in how the parameters
            % are defined
            if contains(fnames{i},'alpha') && sum(y < 0) % check zero crossing of geometrical angle (required for log fit)
                b = 180;
            else
                b = 0;
            end
    
            if contains(fnames{i},'root_offset')
                a = -1;
            else
                a = 1;
            end
    
            y = a*y + b;
    
            normalize = 'on';
            if contains(fnames{i},'width')
                normalize = 'off';
            end
            if strncmp(fit_type,'exp',3) && (numel(x) > 2)
                % a < 0, b < 0
                [cfit{1},gof{1}] = fit(x(:),y(:),'a*exp(b*x)+c','StartPoint',[-mean(y),-1/mean(y),mean(y)],...
                    'Upper',[0,0,Inf],'Lower',[-mean(y),-Inf,0],'Normalize',normalize,'MaxIter',1e3,'TolFun',1e-9);
    
                % a < 0, b > 0
                [cfit{2},gof{2}] = fit(x(:),y(:),'a*exp(b*x)+c','StartPoint',[-mean(y),1/mean(y),mean(y)],...
                    'Upper',[0,1/mean(y),Inf],'Lower',[-Inf,0,0],'Normalize',normalize,'MaxIter',1e3,'TolFun',1e-9);
    
                % a > 0, b < 0
                [cfit{3},gof{3}] = fit(x(:),y(:),'a*exp(b*x)+c','StartPoint',[mean(y),-1/mean(y),mean(y)],...
                    'Upper',[mean(y),Inf,Inf],'Lower',[0,-Inf,0],'Normalize',normalize,'MaxIter',1e3,'TolFun',1e-9);
    
                % a > 0, b > 0
                [cfit{4},gof{4}] = fit(x(:),y(:),'a*exp(b*x)+c','StartPoint',[mean(y),1/mean(y),mean(y)],...
                    'Upper',[Inf,1/mean(y),Inf],'Lower',[0,0,0],'Normalize',normalize,'MaxIter',1e3,'TolFun',1e-9);
                
                [~,idx] = min([cell2mat(gof).sse]);
                cfit = cfit{idx};
    
            elseif strncmp(fit_type,'exp',3) && (numel(x) < 3)
                cfit = fit(x(:),y(:),'exp1');
            elseif strncmp(fit_type,'lin',3)
                cfit = fit(x(:),y(:),'poly1');
            elseif strcmp(fit_type,'poly3')
                cfit = fit(x(:),y(:),'poly3');
                    
            else
                warning('Insufficient data points for exponential fit');
                cfit = fit(x(:),y(:),'poly1');
            end
    
            y = a*y - b;
    
            Prop_fit.(fnames{i}) = @(x) ((a*feval(cfit,x)) - b);
            
            clear cfit;
        end

        if options.plot    
            name = strrep(name,'_',' ');
            name = strrep(name,'alpha','Angle');
            t = uitab(tabgp,'Title',name,'BackgroundColor','w');
            ax = nexttile(tiledlayout('flow',Parent=t));
            hold(ax,'on');
            plot(x, y,'.','Color','k','Parent',ax);
            xx = linspace(0,12,1201);
            plot(xx, feval(Prop_fit.(fnames{i}),xx),'-','Color','r','Parent',ax);
            xlabel('Location [mm]');
            xlim([0 12]);
            yy = ylim;
            if diff(yy) < 10
                ylim(mean(yy)+[-10,10]);
            end
            if startsWith(fnames{i},'r_')
                ylim([0 1]);
            end
            ylabel(name);
            % DisplayUtils.naxis(ax);
            drawnow;
        end
    end
end

