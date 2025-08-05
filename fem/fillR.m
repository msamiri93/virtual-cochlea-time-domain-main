function [partsR,cnt] = fillR(rdof,cdof,m,partsR,cnt)
% generic rectangular matrix fill. degrees of freedom are given to the function as rdof and cdof 

    % row indices
    row = kron(rdof(:),ones(1,length(cdof)));
    row = reshape(row,[],1);    % notice use of reshape

    % column indices
    col = kron(cdof(:),ones(length(rdof),1));

    % Terms (from a single element matrix)
    m = reshape(m,[],1);

    % Ignore if the whole element is restricted by BCs. Such elements have no
    % reduced degrees of freedom, they are all restricted.
    if ~isempty(col)
    
        % This line updates the counter by adding the number of non-zero
        % elements to the last value of the counetr from the previous pass
        cnt = cnt(end) + (1:length(col));
        
        % Save rows, columns and terms coresponding to the current element
        partsR(cnt,1) = row;
        partsR(cnt,2) = col;
        partsR(cnt,3) = real(m);
        partsR(cnt,4) = imag(m);
    
    end

end