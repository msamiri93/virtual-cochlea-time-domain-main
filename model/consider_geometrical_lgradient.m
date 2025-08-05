function gMP = consider_geometrical_lgradient(MP)

    
    A = model_geometry2021_test(MP, 'a');
    B = model_geometry2021_test(MP, 'b');
    
    
    if isfield(MP,'zz')
        zz = MP.zz;
    else
        zz = (-0.5*MP.length_BM:MP.dZ:0.5*MP.length_BM)';
    end
        
%     if isnumeric(MP.loc)    
        if MP.loc<=6       % values will be referenced to the basal properties at Z = 2mm
            zz = 1e3*(MP.loc-2) + zz;            
        else                % values will be referenced to the apiaal properties at Z = 10mm
            zz = 1e3*(MP.loc-10) + zz;
        end
%     end
                                        
    if abs(MP.opt_lgradient) >= 1       % no gradient, the value at z = MP.opt_lgraident is used for the entire span
        zz = MP.opt_lgradient*zz;       % For example,  lgradient=1, forward gradient; lgradient=-1 means reverse gradient
    else
        zz = median(zz)*ones(size(zz));
    end   
                                 
                                        
    POC = {'k_OHB', 'thick_BMa', 'thick_BMp', 'thick_BMz', 'thick_fluid','rOCmass', 'thick_RL', 'mthick_RL', 'thick_TCz',...
           'width_BM', 'diam_DCb', 'diam_DCp', 'thick_TM', 'thick_TMa',  'thick_TMb', 'thick_TMza',  'thick_TMzb', 'width_TM', ...
           'height_TMSLa', 'thick_PC', 'thick_IPC','thick_OPC', 'width_TC', 'height_TC','DC_root_location','OPC_root_location','height_HB', ...
           'angle_RC', 'diam_OHC', 'thick_BM','DD_alpha','DD_theta','diam_DCpr'};       

    gscale = MP.gscale;
    for ii=1:length(POC)
        prop = POC{1,ii};
        if isfield(A,prop)

            grad = set_longi_gradient(A,B,prop,gscale);            
            
            if strcmpi(gscale(1:3),'LIN')
                gMP.(prop) = MP.(prop) + grad*zz;
            elseif strcmpi(MP.gscale(1:3),'LOG')
                gMP.(prop) = MP.(prop)*exp(grad*zz);
            else
                gMP.(prop) = MP.(prop)*exp(0*zz);
            end
            
            if min(gMP.(prop))<0
                wstring = {'Geometric gradient inappropriate';...
                    ['Resuring in negative value of [' prop ']']; ...
                    'Gradient ignored'};
                wname = 'At consider_geometrical_lgradient.m';
                warndlg(wstring, wname);
                gMP.(prop) = mean(gMP.(prop))*ones(size(zz));
            end
            
        end
    end        
    
end

    