function klogical = index2logical(kidx,dof)

klogical = sparse(dof,1);
klogical(kidx)=1;
klogical = logical(klogical);

end