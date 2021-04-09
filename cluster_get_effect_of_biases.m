%% Script for cluster: get effect of biases on biodiversity trends
% Returns estimates of biodiversity change in randomly-selected, abundance-
% biased, and community-biased sites for five change regimes:
% no change, richness increase, richness decrease, total abundance
% increase, and evenness increase

function cluster_get_effect_of_biases(outputdir, JOB_NAME, JOB_ID, change_rate, nsamples, c, g, Stot, mu, var, tend, timelag, start_coordinates)

% Parse input arguments from strings to numbers
change_rate = str2num(change_rate); % relative change rate
nsamples    = str2num(nsamples);
c           = str2num(c);
g           = str2num(g);
Stot        = str2num(Stot);
mu          = str2num(mu);
var         = str2num(var);
tend        = str2num(tend);
timelag     = str2num(timelag);

% Print Toolboxes
ver()

% Get current date
date_current = datestr(now, 'dd/mm/yy-HH:MM');

doAR    = 0;


%% Generate results for different change types

out_none = get_effect_of_biases(c, g, Stot, tend,...
    nsamples, 'none', 0, timelag, doAR, 'mu', mu, 'var', var, start_coordinates);
out_Sinc = get_effect_of_biases(c, g, Stot, tend,...
    nsamples, 'S', change_rate, timelag, doAR, 'mu', mu, 'var', var, start_coordinates);
out_Sdec = get_effect_of_biases(c, g, Stot, tend,...
    nsamples, 'S', -change_rate, timelag, doAR, 'mu', mu, 'var', var, start_coordinates);
out_Ninc = get_effect_of_biases(c, g, Stot, tend,...
    nsamples, 'N', change_rate, timelag, doAR, 'mu', mu, 'var', var, start_coordinates);
out_einc = get_effect_of_biases(c, g, Stot, tend,...
    nsamples, 'evenness', change_rate, timelag, doAR, 'mu', mu, 'var', var, start_coordinates);

fprintf('\n\nFinished successfully.')


%% Generate logfile and save results to workspace

% Save to workspace
workspacefilename = sprintf('%s/workspace_%s_%s', outputdir, JOB_NAME, JOB_ID);
save(workspacefilename)

% Generate logfile
logfilename = sprintf('%s/%s_%s.out', outputdir, JOB_NAME, JOB_ID);
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

