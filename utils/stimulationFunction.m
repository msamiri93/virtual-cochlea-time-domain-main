function y = stimulationFunction(amplitude,time,options)
    arguments
        amplitude       (1,:) {mustBeFloat}
        time
        options.type          {mustBeMember(options.type,{'impulse','tone','random'})} = 'random';
        options.freq    (1,:) {mustBeFloat}
        options.phi     (1,:) {mustBeFloat}
        options.seed (1,1) {mustBeNonnegative,mustBeInteger} = 42;
    end

    rng(options.seed);

    switch options.type
        case 'tone'
            nf = numel(options.freq);
            if ~isfield(options,'phi')
                options.phi = 2*pi*rand(1,nf);
            end
            y = toneStimulationFunction(amplitude, options.freq, options.phi, time);
        case 'random'
            y = randomStimulationFunction(amplitude,options.freq,time);
        case 'impulse'
            y = impulseStimulationFunction(amplitude,time);
    end

end

function y = impulseStimulationFunction(amplitude,t)
    impulse_span = 0.2; % ms
    if t < impulse_span % ms
        h = 0.5*(1 - cos(2*pi*t/impulse_span));
    else
        h = 0;
    end
    y = h*amplitude; % prescribed pressure applied       
end

function y = toneStimulationFunction(amplitude,freq, phi, t)
    nf = numel(freq);
    na = numel(amplitude);
    if isempty(phi)
        phi = 2*pi*rand(1,nf);
    end
    if eq(na,1)
        amplitude = amplitude/sqrt(na);
    end
    y = rampFunction(freq,t)*amplitude.*real(exp(1i*(2*pi*freq*t + phi)));
end

function y = randomStimulationFunction(amplitude,freq,t)
% Generate random frequencies and phases
    lower_bound = 0.5;
    upper_bound = 30;
    num_freq = 121;
    a = log10(lower_bound);
    b = log10(upper_bound);
    ff = 10.^(a + (b - a).*rand(num_freq,1));
    phi = 2*pi*rand(num_freq,1);

    amplitude = amplitude/sqrt(numel(ff));
    y = rampFunction(freq,t)*amplitude*sum(real(exp(1i*(2*pi*ff*t + phi))));
end

function h = rampFunction(freq,t)
    % half-sin initial
    raise_time = 5/freq;
    if t < raise_time
        h = 0.5*(1 - cos(2*pi*t/(2*raise_time)));
    % elseif t > decay_time
    %     h = 0.5*(1 + cos(2*pi*time/(2*MP.period)));
    else
        h = 1;
    end
end