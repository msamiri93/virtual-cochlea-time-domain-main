function Prop_fit = eval_damping_coef(Prop_fit)

    xi = linspace(0.5,11.5,21);
    % r = zeros(size(xi));
    % for i = 1:numel(xi)
    %     k = virtualOoCStatic10(xi(i),lambdaF=-1,update_coef=Prop_fit); 
    %     r(i) = k.OCC_BMC;
    % end
    % 
    % Prop_fit.OOC_stiffness = @(x) (interp1(xi,r,x,"pchip","extrap"));
    % Prop_fit.alpha_c = @(x) (0.5./interp1(xi,r,x,"pchip","extrap"));

    Prop_fit.alpha_c = @(x) (0.05./interp1(xi,2*pi*loc2freq(xi),x,"pchip","extrap"));
    % Prop_fit.alpha_c = @(x) (0.15./interp1(xi,2*pi*loc2freq(xi),x,"pchip","extrap"));
    % Prop_fit.beta_c = @(x) (4*0.01.*interp1(xi,2*pi*loc2freq(xi),x,"pchip","extrap"));
    Prop_fit.beta_c = @(x) (9.5*ones(size(x)));
    % Prop_fit.beta_c = @(x) (0.01*ones(size(x)));
end