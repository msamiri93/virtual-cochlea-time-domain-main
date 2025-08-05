function run_0509

% ff = logspace(log10(3),log10(30),32); % 16/oct preferred
ff = logspace(log10(0.5),log10(25),8);

sensitivity = [0,1];
SPL = 40;
n_rows = 1;

parameters = combinations(SPL,ff,sensitivity);
np = size(parameters,1);

cores = 32;

numWorkers = min([np, feature('numcores')]);
delete(gcp('nocreate'));
parpool('local', numWorkers);

% proj = openProject('Virtualcochleatimedomain.prj');

% kA = logspace(-0.4,0.4,5);
% kA = 1;
% for ii = 1:numel(kA)
% 
%     param.kA = kA(ii);
    % profile on 
    spmd

        % Create a distributed array initialized to cell (or whatever type you need)
        D = codistributed.cell(np, 1);  % 1D column array distributed among workers

        maxNumCompThreads(floor(cores/numWorkers));  

        % Get local index range for this worker
        localRange = globalIndices(D, 1);  % index range this worker is responsible for
        % Populate local part
        localPart = cell(size(getLocalPart(D)));  % same size as local chunk
        
        for k = 1:length(localRange)
            
            % globalIdx = localRange(k);
            % fprintf('Worker %d: %d threads available | Running for : sensitive cochlea: %d, frequency %2.2f kHz, SPL: %2.2f dB\n', ...
            %     spmdIndex, maxNumCompThreads,parameters(globalIdx,:).sensitivity,parameters(globalIdx,:).ff,parameters(globalIdx,:).SPL);
            % 
            % localPart{k} = virtualCochleaDynamic(freq = parameters(globalIdx,:).ff, fluid_Corti = 1, ...
            %     sensitive_cochlea = parameters(globalIdx,:).sensitivity, nonlin = 1, ...
            %     SPL = parameters(globalIdx,:).SPL, cycles = 20, save = true);

            globalIdx = localRange(k);
            fprintf('Worker %d: %d threads available | Running for : sensitive cochlea: %d, frequency %2.2f kHz, SPL: %2.2f dB\n', ...
                spmdIndex, maxNumCompThreads,parameters(globalIdx,:).sensitivity,parameters(globalIdx,:).ff,parameters(globalIdx,:).SPL);
    
            localPart{k} = virtualCochleaDynamic(freq = parameters(globalIdx,:).ff, fluid_Corti = 0, ...
                sensitive_cochlea = parameters(globalIdx,:).sensitivity, nonlin = 1, ...
                SPL = parameters(globalIdx,:).SPL, cycles = 100, save = true, stim_type = 'tone');
        end

        % Assign local result back to distributed array
        D = codistributed.build(localPart, getCodistributor(D));
        
    end

    tmp = gather(D);
    R = cat(1,tmp{:}); 
    
    ct = clock;
    date_str = [num2str(ct(2),'%02.2d'),num2str(ct(3),'%02.2d'),num2str(rem(ct(1),2000),'%02.2d')];
    time_str = [num2str(ct(4),'%02.2d'),num2str(ct(5),'%02.2d')];
    save(['.\houtput\result_',date_str,'_',time_str,'.mat'],'R','-mat','-v7.3');
% end

end

function pulled = random_pull_unique(start_val, end_val)
    % Generate the full range
    if start_val <= end_val
        vals = start_val:end_val;
    else
        vals = start_val:-1:end_val;
    end

    % Shuffle randomly
    shuffled = vals(randperm(length(vals)));
    
    pulled = shuffled;  % Return pulled list if needed
end