%% Vary the input parameter ("parameter") in a specified range ("v_vec")

function [out, grain2mean] = get_parameter_variation(Stot, c, g, tend, nsamples, change, change_rate, v_vec, parameter, doAR, varargin)

assert(isequal(size(Stot),[1,1]), 'First input argument must be double defining total richness')
assert(ismember(change, {'none', 'N', 'S', 'evenness'}), 'Undefined change type.')

% if specified, mean of lognormal SAD distribution
mu = 3; % warning: might not be up to date with default
if any(strcmp(varargin, 'mu'))
    ind = find(strcmp(varargin, 'mu'));
    mu = varargin{ind+1};
end

% if specified, variance of lognormal SAD distribution
var = 1.5;  % warning: might not be up to date with default
if any(strcmp(varargin, 'var'))
    ind = find(strcmp(varargin, 'var'));
    var = varargin{ind+1};
end

% timelag (only every "timelag"st point is used for final timeseries)
timelag = 0;
if any(strcmp(varargin, 'timelag'))
    ind = find(strcmp(varargin, 'timelag'));
    timelag = varargin{ind+1};
end

% spatial distribution of start coordinates
start_coordinates = 'random';
if any(strcmp(varargin, 'aggregated'))
    start_coordinates = 'aggregated';
end

% get spatial heterogeneity? (difference of cell to mean in landscape)
get_heterogeneity = 'off';
if any(strcmp(varargin, 'get_heterogeneity'))
    get_heterogeneity = 'on';
end

% Default: no truncating. Alternatively: truncate pt first time points of timeseries
truncate   = 'off';
pt = 0;

v_num = length(v_vec);

out.var = v_vec;

%% Start variation

tend_input = tend;

for i = 1:v_num
    %% Vary the parameter
    
    fprintf('\nvary %s, v = %d/%d', parameter, i, v_num)
    
    if strcmp(parameter, 'g')
        g = v_vec(i);
    elseif strcmp(parameter, 'c')
        c = v_vec(i);   
    elseif strcmp(parameter, 'cg_factor')
        c = round(c*v_vec(i));   
        g = round(g*v_vec(i));
    elseif strcmp(parameter, 'Stot')
        Stot = v_vec(i); 
    elseif strcmp(parameter, 'mu')
        mu = v_vec(i);   
    elseif strcmp(parameter, 'var')
        var = v_vec(i);       
    elseif strcmp(parameter, 'tend')
        tend = v_vec(i);  
    elseif strcmp(parameter, 'timelag')
        timelag = v_vec(i);    
    elseif strcmp(parameter, 'pct_trunc')
        pct_trunc = v_vec(i);    
        pt = round((pct_trunc/100)*tend_input);
        tend = tend_input + pt; % make timeseries longer, such that after truncating 
        % the intitially provided number of points remain
    elseif strcmp(parameter, 'change_rate')
        change_rate = v_vec(i);
    else
        error('Parameter not yet defined for variation.')
    end
    
    %% Get results for the different biases

    [timeseries_random, grain2mean_random] = get_diversity_timeseries(Stot, c, g, tend,...
        nsamples, change, change_rate, 'timelag', timelag,...
        'mu', mu, 'var', var, start_coordinates, 'truncate', pt);
    if doAR == 1 % if auto-regression is on, include it in errors of regression
        [trends_random]        = get_trends_AR(timeseries_random);
    else
        [trends_random]        = get_trends(timeseries_random);
    end
    samples = trends_random;
    fields  = fieldnames(samples);
    for f = 1 : length(fields)
        field = eval(sprintf('samples.%s', fields{f}));
        if isequal(size(field), [nsamples 1])
            eval( sprintf('out.%s_random(i,:) = samples.%s;', fields{f}, fields{f}) )
        elseif isequal(size(field), [1 1])
            eval( sprintf('out.%s_random(i) = samples.%s;', fields{f}, fields{f}) )
        end
    end
    
    [timeseries_comm_biased, grain2mean_comm_biased] = get_diversity_timeseries(Stot, c, g, tend,...
        nsamples, change, change_rate, 'sampling', 'comm-biased',...
        change, change_rate, 'timelag', timelag,...
        'mu', mu, 'var', var, start_coordinates, 'truncate', pt);
    if doAR == 1
        [trends_comm_biased]     = get_trends_AR(timeseries_comm_biased);
    else
        [trends_comm_biased]     = get_trends(timeseries_comm_biased);
    end
    samples = trends_comm_biased;
    fields  = fieldnames(samples);
    for f = 1 : length(fields)
        field = eval(sprintf('samples.%s', fields{f}));
        if isequal(size(field), [nsamples 1])
            eval( sprintf('out.%s_comm_biased(i,:) = samples.%s;', fields{f}, fields{f}) )
        elseif isequal(size(field), [1 1])
            eval( sprintf('out.%s_comm_biased(i) = samples.%s;', fields{f}, fields{f}) )
        end
    end
    
    [timeseries_rich_biased, grain2mean_rich_biased] = get_diversity_timeseries(Stot, c, g, tend,...
        nsamples, change, change_rate, 'sampling', 'rich-biased',...
        change, change_rate, 'timelag', timelag,...
        'mu', mu, 'var', var, start_coordinates, 'truncate', pt);
    if doAR == 1
        [trends_rich_biased]     = get_trends_AR(timeseries_rich_biased);
    else
        [trends_rich_biased]     = get_trends(timeseries_rich_biased);
    end
    samples = trends_rich_biased;
    fields  = fieldnames(samples);
    for f = 1 : length(fields)
        field = eval(sprintf('samples.%s', fields{f}));
        if isequal(size(field), [nsamples 1])
            eval( sprintf('out.%s_rich_biased(i,:) = samples.%s;', fields{f}, fields{f}) )
        elseif isequal(size(field), [1 1])
            eval( sprintf('out.%s_rich_biased(i) = samples.%s;', fields{f}, fields{f}) )
        end
    end
    
    %% Get difference of selected grain to mean grain in landscape
    
    grain2mean.N_random(i,:) = grain2mean_random.N;
    grain2mean.N_st_random(i,:) = grain2mean_random.N_st;
    grain2mean.S_random(i,:) = grain2mean_random.S;
    grain2mean.S_st_random(i,:) = grain2mean_random.S_st;
    
    grain2mean.N_comm_biased(i,:) = grain2mean_comm_biased.N;
    grain2mean.N_st_comm_biased(i,:) = grain2mean_comm_biased.N_st;
    grain2mean.S_comm_biased(i,:) = grain2mean_comm_biased.S;
    grain2mean.S_st_comm_biased(i,:) = grain2mean_comm_biased.S_st;
    
    grain2mean.N_rich_biased(i,:) = grain2mean_rich_biased.N;
    grain2mean.N_st_rich_biased(i,:) = grain2mean_rich_biased.N_st;
    grain2mean.S_rich_biased(i,:) = grain2mean_rich_biased.S;
    grain2mean.S_st_rich_biased(i,:) = grain2mean_rich_biased.S_st;
   
end

end