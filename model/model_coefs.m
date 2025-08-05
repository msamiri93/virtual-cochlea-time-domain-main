function OP = model_coefs
                               %[   B   ,  A    ];
    % cf = 'k_OHB'; P.(cf) =     [      0.4,      0.6];
    cf = 'k_OHB'; P.(cf) =     [      0.7,      0.7];

    %cf = 'thick_BM'; P.(cf) =  [      2.05,  3.25, 2.00, 0.7, 0.65]; %0514 results
    cf = 'thick_BM'; P.(cf) =  [    0.7, 1.7, 1.6,0.75];
    cf = 'thick_RL'; P.(cf) =  [      1,      0.8];
    cf = 'diam_DCb'; P.(cf) =  [      1,      1.2];
    cf = 'diam_OHC'; P.(cf) =  [      0.9,      0.6];
    cf = 'thick_TMa'; P.(cf) = [      1,      1];
    cf = 'thick_TMb'; P.(cf) = [      1,      1];

    cf = 'DD_theta'; P.(cf) = [      1.5,      2];

    cf = 'Y_BM'; P.(cf) =      [      1,      1];
    cf = 'Y_OHC'; P.(cf) =     [      1,      1];

    cf = 'thick_IPC'; P.(cf) = [      1.9,      1.1];
    cf = 'thick_OPC'; P.(cf) = [      0.8,      0.4];

    cf = 'Y_OPC'; P.(cf) =     [    0.5*0.02,    0.5*0.01];
    cf = 'Y_IPC'; P.(cf) =     [    1*0.1,    1*0.1];

    cf = 'Y_OPCb'; P.(cf) =    P.('Y_OPC');
    cf = 'Y_IPCb'; P.(cf) =    P.('Y_IPC');

    cf = 'Y_DCb'; P.(cf) =     [  6.5*0.5,    6.5*2];
    % cf = 'Y_DCp'; P.(cf) =     [   0.02,   0.04];
    cf = 'Y_DCp'; P.(cf) =     [   0.005,   0.005];
    
    cf = 'Y_RLx'; P.(cf) =     [      0.5*1,      0.5*1];
    cf = 'Y_RLp'; P.(cf) =     [      0.5*1,      0.5*1];
    
    cf = 'Y_TMx'; P.(cf) =     [      0.8,      1.2];
    % cf = 'Y_TMx'; P.(cf) =     [      40,      50]; % prodanovic

    cf = 'Y_TMz'; P.(cf) =     [    1,   0.5];
    cf = 'Y_YMz'; P.(cf) =     [    1,   0.25];
    cf = 'Y_BMz'; P.(cf) =     [    1,   0.5];
    cf = 'width_BMP'; P.(cf) = [    0.90,   1.05];
    cf = 'width_BMA'; P.(cf) = [    0.90,   1.05];
    % cf = 'tilt_OHC'; P.(cf) =  [    1.0,     1.0];
    % cf = 'r_IPC'; P.(cf) =     1.3;

    cf = 'alpha_c'; P.(cf) =    [   0.85,  0.75,  0.7];
    % cf = 'beta_c'; P.(cf) =     [   0.48,  0.59,  1.14];
    % cf = 'beta_c'; P.(cf) =     [   0.38,  0.36,  0.32];
    cf = 'beta_c'; P.(cf) =     [   1,  1,  1];

    O_names = fieldnames(P);
    
    for i = 1:numel(O_names)
        OP(i).name = O_names{i};
        OP(i).value = P.(O_names{i});
        if numel(OP(i).value) == 2
            OP(i).location = [2,10];
        elseif numel(OP(i).value) == 3
            OP(i).location = [2,4,10];                
        elseif numel(OP(i).value) == 4
            OP(i).location = [0,2,10,12];
        else
            warning(['location for parameter ',O_names{i}, ' was not defined. \n']);
        end
        
    end
    
    save('./assets/optimization_params_0807.mat','OP');

end