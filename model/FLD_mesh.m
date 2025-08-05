function FLD = FLD_mesh(MP,H,FLD)

    H_05 = round(H/2);
    H_OCC = 25;
    L = MP.length_BM;
    h0 = 10;

    if MP.vMC == 1
        load('./assets/vMC_4mm.mat','p','t');
        Nd.x = p; El.node = t;
        El.type = 3;
        Lo = 1000;
        Ho = 80;
        Xo = 2500;
        Ls = 1000;
        Hs = 200;
        Xs = 2500;
        FLD.Lo = Lo;
        FLD.Xo = Xo;
        FLD.Ho = Ho;
        FLD.Ls = Ls;
        FLD.Xs = Xs;
        FLD.Hs = Hs;
    else
        % using quadriateral mesh
        fname = "./assets/SF_12mm_coarse";

        % pyenv('Version','C:\Users\Mohammad Shokrian\AppData\Local\Programs\Python\Python39\python.exe');        
        % pyrunfile(sprintf("./model/gen_gmsh.py " + ...
        %     "--file_name %s --depth_scaling %f --heli_scaling %f " + ...
        %     "--dh %f",fname,1.15,0.975,14));
        % % % % % using quadriateral mesh
        % % % % % fname = "./hinput/SF_12mm_coarse";
        % % % % % pyrunfile(sprintf("./lib/gen_gmsh.py " + ...
        % % % % %     "--file_name %s --depth_scaling %f --heli_scaling %f " + ...
        % % % % %     "--dh %f",fname,1,1,10));        
        t = load(fname);
        t = remove_unused_mesh(t);
        % load('./hinput/mesh_12mm_cochlea_quad_3.mat','msh');
        % Nd.x = msh.POS(:,[1,2]);    
        % El.node = msh.QUADS(:,1:4);
        Nd.x = t.nodes(:,[1,2]);
        El.node = t.elements;
        Nd.bnd = t.boundaries;
        El.type = 4;     
    end

    El.N = size(El.node,1);
    Nd.N = size(Nd.x,1);
    FLD.Nd = Nd;
    FLD.El = El;

    FLD.H = H;
    FLD.h0 = h0;
    FLD.HOC = H_OCC;
    FLD.L = L;

end

function t = remove_unused_mesh(t)

    isValid = (ismember(1:size(t.nodes,1), t.elements(:)));
    mapping = zeros(size(t.nodes,1),1);
    mapping(isValid) = 1:sum(isValid);
    t.elements = mapping(t.elements);
    t.nodes = t.nodes(isValid, :);
    bnames = fieldnames(t.boundaries);
    for i = 1:numel(bnames)
        t.boundaries.(bnames{i}) = mapping(t.boundaries.(bnames{i})); 
    end
end