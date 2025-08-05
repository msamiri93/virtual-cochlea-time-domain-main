function [C, El] = createC2(Nd,El,MP)

global coef;

NDOF = 6;

% Initiate matrix that hold indices of rows and columns and terms of the
% global matrix. First column holds the indices of rows, the second column
% indices of columns, and the third one the terms. 144 comes from the safe
% assumption that all elements have 12 DoF (beam elements). Zero elements
% will be removed at the end. This 3-column matrx is the input for the
% sparse function.
partsC = zeros(144*El.N,3);

% alpha_c = MP.alpha_c*exp(MP.alpha_cgrad*Nd.Z(El.Nd1));
% beta_c = MP.beta_c*exp(MP.beta_cgrad*Nd.Z(El.Nd1));
a = 1e-3; % 1e-3 coef is to translate [um] to [mm] for interpolation
alpha_c = MP.Prop_fit.alpha_c(Nd.Z(El.Nd1)*a + MP.loc);
beta_c = MP.Prop_fit.beta_c(Nd.Z(El.Nd1)*a + MP.loc);

cnt = 0;  % counter variable

Slave_node = false(Nd.N);
Slave_node(Nd.Master_Slave(:,2)) = true;

for i = 1:El.N
    nd1 = El.Nd1(i);
    nd2 = El.Nd2(i);
    if El.type(i) < 1        
        % for the beam elements
%         Ke = beam_eK(El,i);
%         Me = beam_eM(El, i);
        Ke = El.mat(i).k;
        Me = El.mat(i).m;
        if strcmpi(El.name(i),'OHB'),
            alpha = alpha_c(i);
            beta = beta_c(i);
        elseif strcmpi(El.name(i),'OHC'),            
%             c0 = 300; % 300 nN-s/m
            c0 = coef.OHC_damping*300; % 300 nN-s/m
            L0 = 20;  % for 20 um-long OHC
            kOHC = El.YM(i)*El.A(i)/El.oL(i);
            Le = El.oL(i);
            alpha = (c0/kOHC)*(Le/L0);
            beta = 0;
        elseif  strcmpi(El.name(i),'ANK'),
            % Viscous friciton btw TM and RL is lumped to the ANK element
            kANK = El.YM(i)*El.A(i)/El.oL(i);
            eRLx = find(strcmpi(El.name,'RLx1'),1);
            Area = 2.0*El.oL(eRLx)*MP.dZ ;
            mu = 0.72; % 0.72*10^-3  Ns/m^2
            alpha = (mu*Area/MP.height_HB)/kANK;
            beta = 0;
        elseif (isfield(coef,'BM_damping_coef') && strncmpi(El.name(i),'BM',2))
            alpha = 0;   % Proportional damping.
            beta = coef.BM_damping_coef*beta_c(i);
        elseif (isfield(coef,'RL_damping_coef') && strncmpi(El.name(i),'RL',2))
            alpha = 0;   % Proportional damping.
            beta = coef.RL_damping_coef*beta_c(i);
        elseif (isfield(coef,'TM_damping_coef') && strncmpi(El.name(i),'TM',2))
            alpha = 0;   % Proportional damping.
            beta = coef.TM_damping_coef*beta_c(i);            
        elseif (isfield(coef,'DC_damping_coef') && strncmpi(El.name(i),'DC',2))
            alpha = 0;   % Proportional damping.
            beta = coef.DC_damping_coef*beta_c(i); 
        elseif (isfield(coef,'PC_damping_coef') && ~isempty(regexp(El.name(i),regexptranslate('wildcard','*PC'), 'once')))
            alpha = 0;   % Proportional damping.
            beta = coef.PC_damping_coef*beta_c(i);            
        else
            alpha = alpha_c(i);   % Proportional damping.
            beta = beta_c(i);
        end
        c = alpha*Ke + beta*Me;        
        ndof1 = (nd1-1)*NDOF + (1:6);
        ndof2 = (nd2-1)*NDOF + (1:6);               
    else
        % for the link elements
        c = link_eC(Nd,El,MP,i);        
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
    [partsC,cnt] = fillC(Nd,ndof1,ndof2,c,partsC,cnt);
    El.mat(i).c = c;
end

% Remove excess zero elements resulting from 144 DoF assumption.
zr = partsC(:,1)==0;
partsC(zr,:) = [];

% Use sparse function to create the matrix
C = sparse(partsC(:,1),partsC(:,2),partsC(:,3),Nd.rdof,Nd.rdof);

end

function [partsC,cnt] = fillC(Nd,ndof1,ndof2,c,partsC,cnt)

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
row = reshape(row,[],1);

% Terms (from a single element matrix)
c = c(tmp2,tmp2);
c = reshape(c,[],1);

% Ignore if the whole element is restricted by BCs. Such elements have no
% reduced degrees of freedom, they are all restricted.
if ~isempty(col)
    
    % This line updates the counter by adding the number of non-zero
    % elements to the last value of the counetr from the previous pass
    cnt = cnt(end) + (1:length(col));
    
    % Save rows, columns and terms coresponding to the current element
    partsC(cnt,1) = row;
    partsC(cnt,2) = col;
    partsC(cnt,3) = c;
    
end

end

function eC = link_eC(Nd, El, MP, ie)

if MP.opt_lgradient==1,
    dZ = Nd.Z(El.Nd1(ie));
else
    dZ = 0;
end

a = 1e-3; % 1e-3 coef is to translate [um] to [mm] for interpolation
alpha = MP.Prop_fit.alpha_c(Nd.Z(El.Nd1(ie))*a + MP.loc);

if strcmp(El.name(ie),'ANK'),
    mu = 0.7; % 0.7 pN.ms/um^2
    wid = 10; % strip width of one column of 3-pack OHBs
    len = El.oL(find(strcmp(El.name,'RLx1'),1))*2.0; % strip witdh of one column of 3-pack OHBs
    
    if strcmpi(El.name(ie-1),'xHB'),
        hgt = El.oL(ie-1);
    else
        hgt = MP.height_HB;
        fprintf(1,'\nWarning:: make_C2: failed to find element xHB \n');
    end
    
    damp = mu*wid*len/hgt;
    
else
    damp = alpha*El.YM(ie)*El.A(ie)/El.oL(ie);
end

eC = zeros(6);
eC(1,1) = damp;
eC(4,4) = damp;
eC(1,4) = -1*damp;
eC(4,1) = -1*damp;

R = rotation_matrix(El, ie);
eC = R'*eC*R;

end

function R = rotation_matrix(El, ie)

if abs(El.dir(ie,3)-1)<0.1,
    ez = [1, 0, 0];
    % note that most elements alligned either radial(x) or longitudianl(z) direction
else
    ez = [0, 0, 1];
end

xv = El.dir(ie,:);

yv = cross(ez,xv);
yv = yv/norm(yv);
zv = cross(xv,yv);
zv = zv/norm(zv);

T33 = [xv;yv;zv];
O33 = zeros(3);

R = [T33, O33; O33, T33];

end