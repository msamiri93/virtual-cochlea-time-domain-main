function Prop = apply_model_coefs(Prop)

    P_names = {Prop.name};

    for i = 1:numel(Prop)
        Prop(i).location = [2,10];
        Prop(i).value = [1,1];
    end

    load('./assets/optimization_params.mat','OP');
    O_names = {OP.name};

    for i = 1:numel(OP)

        if strcmp('Y_OCC',O_names{i}) % adjusting OoC YM longitudinal gradient
            idx = find(strncmp(P_names,'Y_',2) & ~contains(P_names,'z') ...
                & ~contains(P_names,'BM') & ~contains(P_names,'HB') ...
                & ~contains(P_names,'ANK') & ~contains(P_names,'OHC') ...
                & ~contains(P_names,'PC'));
        elseif strcmp('Y_YMz',O_names{i}) % adjusting YM longitudinal stiffness
            idx = find(strncmp(P_names,'Y_',2) & contains(P_names,'z') ...
                & ~contains(P_names,'BM') & ~contains(P_names,'TM'));
        elseif strcmp('Y_BMz',O_names{i}) % adjusting YM longitudinal stiffness
            idx = find(strncmp(P_names,'Y_',2) & contains(P_names,'z') & contains(P_names,'BM'));                
        elseif strcmp('thick_OCC',O_names{i}) % adjusting OoC elements' thickness gradient
            idx = find(contains(P_names,'thick') & contains(P_names,'diam') ...
                & ~contains(P_names,'z') & ~contains(P_names,'BM') ...
                & ~contains(P_names,'OHC'));           
        else
            idx = find(contains(P_names,O_names{i}) & ~contains(P_names,'z'));
        end

        for j = 1:numel(idx)
            if isempty(OP(i).location)
                Prop(idx(j)).value = OP(i).value;
            else
                Prop(idx(j)).value = OP(i).value;
                Prop(idx(j)).location = OP(i).location;
            end
        end
    end

end