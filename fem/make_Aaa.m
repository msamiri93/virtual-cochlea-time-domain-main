function [Aaa,Iaa] = make_Aaa(MP,HB)   
    dof = MP.dof(4).n;
    Aaa = spdiags(reshape([HB.kA].*[HB.kES],dof,1),0,dof,dof);
    % an i*w*I, where I is the identity matrix is needed at frequency
    % dependent assembly step
    Iaa = speye(dof);
end