%% Script for cluster: get effect of biases on biodiversity trends
% Returns estimates of biodiversity change in randomly-selected, abundance-
% biased, and community-biased sites for three change regimes
% no change, richness increase, and richness decrease, across gradients in
% either 

function cluster_get_scale_dependency(outputdir, JOB_NAME, JOB_ID, change_rate, nsamples, variable, nvar, start_coordinates)

% Parse input arguments from strings to numbers
nsamples      = str2num(nsamples);
change_rate   = str2num(change_rate); % relative change rate
nvar          = str2num(nvar); % number of grain sizes

% Print Toolboxes
ver()

% Get current date
date_current = datestr(now, 'dd/mm/yy-HH:MM');

%% Default parameterization 

% Landscape (c = length of quadratic grid, g = length of quadratic grain)
c = 100;
g = 5;

% Species abundance distribution parameters
Stot = 100; 
mu   = 4;
var  = 1.5;

% Simulation parameters
tend    = 20;
timelag = 0;
doAR    = 0;

%% Generate results

% define variation vector for the chosen variable
switch variable
    case 'g' % vary grain size
        % vary g such that relative g (to c) will be equidistant
        % % v_vec = round(linspace(1, 50, nvar));
        gsquare_vec = linspace(0.1, 1, nvar);
        v_vec = unique(round(sqrt(gsquare_vec*c*c/100)));
    case 'gfull' % vary grain size (wider range of g)
        gsquare_vec = [0.1 0.25 0.5 1 2 5 10 20 50];
        v_vec = round(sqrt(gsquare_vec*c*c/100));
        variable = 'g';
    case 'c' % vary landscape size
        csquare_vec = linspace(50*50, 250*250, nvar);
        v_vec = g*round(sqrt(csquare_vec)/g); % round to next multiple of g
    case 'cg_factor' % vary factor between landscape and grain size
        v_vec = linspace(0.5, 2.5, nvar);    
    case 'Stot' % vary total number of species
        v_vec = round(linspace(50, 150, nvar));
    case 'mu' % vary mean of SAD
        v_vec = linspace(3, 5, nvar); 
    case 'var' % vary variance of SAD
        v_vec = linspace(1, 2, nvar);
    case 'tend' % vary sampling duration
        v_vec = round(linspace(5, 50, nvar));    
    case 'timelag' % vary time lag between sampling events
        v_vec = round(linspace(0, 10, nvar));
    case 'pct_trunc' % vary percentage of time series truncated
        v_vec = 0:5:50;
    case 'change_rate' % vary rate of biodiversity change
        v_vec = linspace(0, 0.0025*10, nvar);    
        
end
        
% No change
if ~strcmp(variable, 'change_rate')
    change = 'none';
    [var_none, grain2mean_none] = get_parameter_variation(Stot, c, g, tend, nsamples,...
        change, 0, v_vec, variable, doAR, 'mu', mu, 'var', var, 'timelag', timelag, start_coordinates);
    fprintf('\nfinished variation of "no change".')
end

% S increase
change = 'S';
[var_Sinc, grain2mean_Sinc] = get_parameter_variation(Stot, c, g, tend, nsamples,...
    change, change_rate, v_vec, variable, doAR, 'mu', mu, 'var', var, 'timelag', timelag, start_coordinates);
fprintf('\nfinished variation of "S increase".')

% S decrease
change = 'S';
[var_Sdec, grain2mean_Sdec] = get_parameter_variation(Stot, c, g, tend, nsamples,...
    change, -change_rate, v_vec, variable, doAR, 'mu', mu, 'var', var, 'timelag', timelag, start_coordinates);
fprintf('\nfinished variation of "S decrease".')

fprintf('\n\n Finished successfully.\n\n')

%% Generate logfile and save results to workspace

% Save to workspace
workspacefilename = sprintf('%s/workspace_%s_%s_%s', outputdir, JOB_NAME, JOB_ID, variable);
save(workspacefilename)

% Generate logfile
logfilename = sprintf('%s/%s_%s_%s.out', outputdir, JOB_NAME, JOB_ID, variable);
fprintf(logfilename)
logfile = fopen(logfilename, 'w');
fprintf(logfile, 'doAR: %d\n', doAR);
fprintf(logfile, 'Stot: %d\n', Stot);
fprintf(logfile, 'c: %d\n', c);
fprintf(logfile, 'g: %d\n', g);
fprintf(logfile, 'timelag: %d\n', timelag);
fprintf(logfile, 'tend: %d\n', tend);
fprintf(logfile, 'change_rate: %1.4f\n', change_rate);
fprintf(logfile, 'nsamples: %d\n', nsamples);
fclose(logfile);

end

