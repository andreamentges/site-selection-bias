%% Returns the slope estimates of richness for three site-selection strategies

function [out, pct_change, pct_dec] = get_effect_of_biases(c, g, Stot, tend, nsamples, change, change_rate, timelag, doAR, varargin)

% checks
assert(doAR==0 | doAR==1, 'doAR must be binary.')
assert(ismember(change, {'none', 'N', 'S', 'evenness'}), 'Undefined change type.')
assert(isequal(size(Stot), [1 1]), 'Specify Stot as third input argument.')

% Default: no truncating (use whole time-series).
% Alternatively: truncate pt percent of timeseries
truncate   = 'off';
pt = 0;
if any(strcmp(varargin, 'truncate'))
    truncate = 'on';
    ind = find(strcmp(varargin, 'truncate'));
    pt = varargin{ind+1}; 
end

if change_rate > 0
    fprintf('\nGet results for %s increase (%d samples)...', change, nsamples)
elseif change_rate < 0
    fprintf('\nGet results for %s decrease (%d samples)...', change, nsamples)
else 
    fprintf('\nGet results for no change (%d samples)...', nsamples)
end

% mean of lognormal SAD distribution
if any(strcmp(varargin, 'mu'))
    ind = find(strcmp(varargin, 'mu'));
    mu = varargin{ind+1};
end

% variance of lognormal SAD distribution
if any(strcmp(varargin, 'var'))
    ind = find(strcmp(varargin, 'var'));
    var = varargin{ind+1};
end

% spatial distribution of start coordinates
start_coordinates = 'random';
if any(strcmp(varargin, 'aggregated'))
    start_coordinates = 'aggregated';
elseif any(strcmp(varargin, 'heterogeneous'))
    start_coordinates = 'heterogeneous';
end


%% Get diversity time-series in a sampling grain
         
% Check whether the variables mu and var exist
if exist('mu', 'var') && exist('var', 'var')
    % get timeseries for SAD with specified mu and var
    [timeseries_random] = get_diversity_timeseries(Stot, c, g, tend, nsamples,...
    change, change_rate, 'sampling', 'random', 'timelag', timelag,...
    'truncate', pt, 'mu', mu, 'var', var, start_coordinates);
    [timeseries_comm_biased] = get_diversity_timeseries(Stot, c, g, tend, nsamples,...
        change, change_rate, 'sampling', 'comm-biased', 'timelag', timelag,...
        'truncate', pt, 'mu', mu, 'var', var, start_coordinates);
    [timeseries_rich_biased] = get_diversity_timeseries(Stot, c, g, tend, nsamples,...
        change, change_rate, 'sampling', 'rich-biased', 'timelag', timelag,...
        'truncate', pt, 'mu', mu, 'var', var, start_coordinates);
    
elseif ~exist('var', 'var') && ~exist('var', 'var')
    % get timeseries from default SAD for each sample
    [timeseries_random] = get_diversity_timeseries(Stot, c, g, tend, nsamples,...
    change, change_rate, 'sampling', 'random', 'timelag', timelag,...
    'truncate', pt, start_coordinates);
    [timeseries_comm_biased] = get_diversity_timeseries(Stot, c, g, tend, nsamples,...
        change, change_rate, 'sampling', 'comm-biased', 'timelag', timelag,...
        'truncate', pt, start_coordinates);
    [timeseries_rich_biased] = get_diversity_timeseries(Stot, c, g, tend, nsamples,...
        change, change_rate, 'sampling', 'rich-biased', 'timelag', timelag,...
        'truncate', pt, start_coordinates);
    
else
    error('Did you specify only one of mu and var for SAD?')
end


%% Get biodiversity trend estimates from the time series

if doAR == 0
    out.random      = get_trends(timeseries_random);
    % out.pop_biased  = get_trends(timeseries_pop_biased);
    out.comm_biased = get_trends(timeseries_comm_biased);
    out.rich_biased = get_trends(timeseries_rich_biased); 
else
    error('AR is deactivated in the get_trends_AR.m')
    out.random      = get_trends_AR(timeseries_random);
    % out.pop_biased  = get_trends_AR(timeseries_pop_biased);
    out.comm_biased = get_trends_AR(timeseries_comm_biased);
    out.rich_biased = get_trends_AR(timeseries_rich_biased); 
end



end

