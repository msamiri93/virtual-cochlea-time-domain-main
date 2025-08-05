function k = perm_scale(MP,x)

    x = x*1e-3 + MP.loc;
    kb = MP.kk; % 2mm
    ka = 1e0*MP.kk; % 10mm

    ab_dist = 8;
    k = zeros(length(x),1);
    if MP.kk ~= 0
        grad = log(kb./ka)/ab_dist;
        k = exp(log(kb) + grad*(x-2));
    end

end