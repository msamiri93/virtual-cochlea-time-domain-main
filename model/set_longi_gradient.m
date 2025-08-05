function grad = set_longi_gradient(A, B, prop, pscale)
    
    aval = A.(prop);
    bval = B.(prop);
    ab_dist = 8000;     % Distance between the reference points at base (z = 2mm) and apex (z = 10 mm)
    
    if ~exist('pscale','var')
        pscale = 'LOG';
    end
    
    if norm(bval)~=0
        if strcmpi(pscale(1:3),'LOG')
            grad = log(aval./bval)/ab_dist;
        else
            grad = (aval-bval)/ab_dist;
        end
    else
        grad = zeros(size(A.(prop)));
    end
    
end % of fucntion set_longi_gradient()