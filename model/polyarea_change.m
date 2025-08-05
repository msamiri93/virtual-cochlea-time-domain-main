function H = polyarea_change(xx)
% shoelace formula for infinitesimal deformation
% H is multiplied by [dx_1; ... ;dx_n;dy_1; ... ;dy_n] 
% where x and y displacements are arranged vertically
% in order for the calculation to be positive, the nodal arrangement should
% be counterclock-wise

    x = xx(:,1);
    y = xx(:,2);

    siz = size(xx);
    dx = x([siz(1),1:siz(1)-1],:) - x([2:siz(1),1],:);
    dy = - y([siz(1),1:siz(1)-1],:) + y([2:siz(1),1],:);

    H = 0.5*[dy;dx];
    H = H.';

end