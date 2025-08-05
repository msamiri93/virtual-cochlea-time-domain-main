function [C,G,I,V0] = mElectrical2_cmod_MSA(MP,OHC,IHC)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % This is the analysis of whole three scalae circuit.This is giving output for dynamic condition.
% % This circuit is for solving the differential equation for the circuit formed in the function "tj_OClongnew_24thapril_real_dynamic2"
% % No input for this function explicitly but it reads the BM vibration from the mat file 'input', the output is time and node voltages
% % First written in April,2013 as a part of the course project BME 404
% % The circuit is a combination Mistrik et al. 2009 paper and NAMLAB circuit

   % Initial condition
   
   nz = round(MP.length_BM/MP.dZ)+1;           % # of sections along longitudinal direction
   nv = 5;
   % MSA: these definition for nodes are not clear in where they sit in the
   % matrix for the equations
    % % Voltage vector, V= [ V1,1; V1,2; V1,4; V1,6; V1,7; V1,8; V1,9];
    % % V1,1= Node of scala media
    % % V1,2= Intracellular space of OHC
    % % V1,4= Extracellular space of OHC
    % % V1,6= Intracellular space of IHC
    % % V1,7= Extracellular space of IHC
    % % V1,8= Node of scala tympani
    % % V1,9= Node of scala vestibuli
    % % %%%%%%%%Read input file
    [StV,SL,OoC] = organ_Corti_parameters_wgrad(MP,nz);

    [G0,I]= make_G0matrix(nv,nz,StV,OHC,IHC,SL,OoC);
    V0 = G0\I; % steady state electrial potential
    [C,G] = make_CGmatrix(MP,nv,nz,V0,StV,OHC,IHC,SL,OoC);

end
