function [K,El] = createK(Nd,El)

NDOF = 6;

% Initiate matrix that hold indices of rows and columns and terms of the
% global matrix. First column holds the indices of rows, the second column
% indices of columns, and the third one the terms. 144 comes from the safe
% assumption that all elements have 12 DoF (beam elements). Zero elements
% will be removed at the end. This 3-column matrx is the input for the
% sparse function.
partsK = zeros(144*El.N,3);
cnt = 0; % counter variable

Slave_node = false(Nd.N);
Slave_node(Nd.Master_Slave(:,2)) = true;

for i = 1:El.N
    
    nd1 = El.Nd1(i);
    nd2 = El.Nd2(i);
    
    % for the beam elements
    if El.type(i) < 1        
        k = beam_eK(El,i);
        ndof1 = (nd1-1)*NDOF + (1:6);
        ndof2 = (nd2-1)*NDOF + (1:6);        
    else
        % for the link elements
        k = link_eK(El,i);
        ndof1 = (nd1-1)*NDOF + (1:3);
        ndof2 = (nd2-1)*NDOF + (1:3);
    end

    if El.hinge(i,1) && Slave_node(nd1)
        nd = Nd.Master_Slave(eq(Nd.Master_Slave(:,2),nd1),1);
        ndof1(1:3) = (nd-1)*NDOF + (1:3);
    elseif El.hinge(i,2) && Slave_node(nd2)
        nd = Nd.Master_Slave(eq(Nd.Master_Slave(:,2),nd2),1);
        ndof2(1:3) = (nd-1)*NDOF + (1:3);
    end

    % This function fills in the 3-column input for the sparse
    % function.
    [partsK,cnt] = fillK(Nd,ndof1,ndof2,k,partsK,cnt);
        
    El.mat(i).k = k;
end

% Remove excess zero elements resulting from 144 DoF assumption.
zr = partsK(:,1)==0;
partsK(zr,:) = [];

% Use sparse function to create the matrix
K = sparse(partsK(:,1),partsK(:,2),partsK(:,3),Nd.rdof,Nd.rdof);

end

function [partsK,cnt] = fillK(Nd,ndof1,ndof2,k,partsK,cnt)

% Function updates 3-column input martix and the counter variable.

% Get reduced DoF for the element
tmp = [Nd.mapF2R(ndof1); Nd.mapF2R(ndof2)];
tmp1 = tmp(tmp~=0); % remove zeros (BCs)
tmp2 = find(tmp);   % logical index varable of non-zero terms

% Kron function serves to repeat certain patterns.
% For example: we have [1; 2; 3] and want to get [1; 1; 2; 2; 3; 3;].
% To obtain [1; 1; 2; 2; 3; 3;] we use kron([1;2;3],ones(2,1)). Paired with
% reshape and transpose, as well as chanign order of inputs, kron(a,b) vs.
% kron(b,a), we can get any repeating pattern.

% column indices
col = kron(tmp1,ones(length(tmp1),1));

% row indices
row = kron(tmp1,ones(1,length(tmp1)));
row = reshape(row,[],1);    % notice use of reshape

% Terms (from a single element matrix)
k = k(tmp2,tmp2);
k = reshape(k,[],1);

% Ignore if the whole element is restricted by BCs. Such elements have no
% reduced degrees of freedom, they are all restricted.
if ~isempty(col)
    
    % This line updates the counter by adding the number of non-zero
    % elements to the last value of the counetr from the previous pass
    cnt = cnt(end) + (1:length(col));
    
    % Save rows, columns and terms coresponding to the current element
    partsK(cnt,1) = row;
    partsK(cnt,2) = col;
    partsK(cnt,3) = k;
    
end

end