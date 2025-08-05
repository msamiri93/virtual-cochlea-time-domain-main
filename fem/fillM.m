function [partsM,cnt] = fillM(NV, nodes, m, partsM, cnt)
% written for complex numbers       
    dof = kron(nodes, ones(NV,1)); % NV is the number of variables at each node
    dof = NV*dof - repmat(fliplr(0:NV-1)',length(nodes),1);
    col = kron(dof,ones(length(dof),1));
    row = repmat(dof,length(dof),1);
    m = m(:);

    cnt = cnt(end) + (1:length(col));

    partsM(cnt,1) = row;
    partsM(cnt,2) = col;
    partsM(cnt,3) = real(m);
    partsM(cnt,4) = imag(m);

end