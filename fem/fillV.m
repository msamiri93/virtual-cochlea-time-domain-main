function B = fillV(NV, nodes, F_e, B)
% fill a vector based on number of variables and nodes numbers
% vector of variables format: [u1;v1;p1;w1;...;un;vn;pn;wn];

    dof = kron(nodes, ones(NV,1)); % NV is the number of variables at each node
    
    dof = NV*dof - repmat(fliplr(0:NV-1)',length(nodes),1);
    
    B(dof,1) = B(dof,1) + F_e;
    
end