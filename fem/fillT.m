function [partsT,cnt] = fillT(rdof,cdof,v,partsT,cnt)

    % fill in the Triplets
    cnt = cnt(end) + (1:numel(v));
    partsT(cnt,1) = rdof(:);
    partsT(cnt,2) = cdof(:);
    partsT(cnt,3) = real(v(:));
    partsT(cnt,4) = imag(v(:));
end