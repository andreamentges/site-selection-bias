%% Main script
% Generates the simulation data and plots results.
%
% Computer code corresponding to the following publication:
% A. Mentges, S. Blowes, D. Hodapp, H. Hillebrand, and J. Chase
% "Effects of site-selection bias on estimates of biodiversity change",
% Conservation Biology, 2020
%
% Author of code: Andrea Mentges, 2019

% Default font sizes
set(groot,'defaultAxesFontSize', 11)
set(groot,'defaultTextFontSize', 12)

%% Default set-up

c        = 100;     % size of grid (c x c cells)
g        = 5;       % size of potential sites (g x g cells)
Stot     = 100;     % size of regional species pool
tend     = 20;      % number of time steps
timelag  = 0;       % number of time steps between sampling events
nsamples = 10000;   % number of samples taken (1 per landscape)

% relative change rate (percent of individuals that change in each time
% step for change regimes)
r_change_rel = 0.0025; 

% mean and standard deviation of species abundance distribution
mu       = 4;
sd      = 1.5;

% initial distribution of individuals (default is '' (for random),
% alternatively use 'aggregated')
start_coordinates = 'random'; 

%% Get results (Figure 3): Effect of biases

% Start simulation
output_folder = cd;
job_name      = 'effect_of_biases';
job_ID        = '001';
cluster_get_effect_of_biases(output_folder, job_name, job_ID, ...
    num2str(r_change_rel), num2str(nsamples), num2str(c), num2str(g),...
    num2str(Stot), num2str(mu), num2str(sd), num2str(tend),...
    num2str(timelag), start_coordinates);
load(sprintf('%s/workspace_%s_%s.mat', output_folder, job_name, job_ID))

%% Plot Figure 3: Effect of biases

% choose aggregated or random spatial distribution
sp_dist = 'aggregated'; % for Figure S8
sp_dist = 'random'; % for Figure 3
        
fc = [0.5 0.5 0.5;
    mycolors('lightorange');
    mycolors('redorange')];

yl_col = [.8 .8 .8];
yl_width = 2;

boxwidth = 0.1;

trotation = 60;

figure('color', 'white', 'position', [-781,263,459,265])
subplot(1,3,1)
ax1 = gca;
subdata = [out_none.random.S_slope out_none.comm_biased.S_slope ...
    out_none.rich_biased.S_slope];
fplot_boxplot(subdata, ax1) 
set(gca, 'XAxisLocation', 'origin', 'Box', 'off')
y1 = ylabel(sprintf('Estimated\nrichness change'));
title1 = title('{\bf a}    No change', 'FontWeight', 'bold', ...
    'HorizontalAlignment', 'left');
t1 = text(ax1, [0.1 0.2 0.3], repmat(ax1.YLim(1)*1.1, 1, 3),...
    {'Random', 'Comm. bias', 'Rich. bias'},...
    'Rotation', trotation, 'HorizontalAlignment', 'right', 'Fontsize', 11);

subplot(1,3,2)
ax2 = gca;
subdata = [out_Sdec.random.S_slope out_Sdec.comm_biased.S_slope ...
    out_Sdec.rich_biased.S_slope];
fplot_boxplot(subdata, ax2) 
set(gca, 'XAxisLocation', 'origin', 'Box', 'off', 'XTickLabels', [],...
    'YTickLabels', [])
title2 = title('{\bf b}   S decrease', 'FontWeight', 'bold', ...
    'HorizontalAlignment', 'left');
t2 = text(ax2, [0.1 0.2 0.3], repmat(ax2.YLim(1)*1.1, 1, 3),...
    {'Random', 'Comm. bias', 'Rich. bias'},...
    'Rotation', trotation, 'HorizontalAlignment', 'right', 'Fontsize', 11);

subplot(1,3,3)
ax3 = gca;
subdata = [out_Sinc.random.S_slope out_Sinc.comm_biased.S_slope ...
    out_Sinc.rich_biased.S_slope];
fplot_boxplot(subdata, ax3) 
set(gca, 'XAxisLocation', 'origin', 'Box', 'off', 'XTickLabels', [], ...
    'YTickLabels', [])
title3 = title('{\bf c}    S increase', 'FontWeight', 'bold', ...
    'HorizontalAlignment', 'left');
t3 = text(ax3, [0.1 0.2 0.3], repmat(ax3.YLim(1)*1.1, 1, 3),...
    {'Random', 'Comm. bias', 'Rich. bias'},...
    'Rotation', trotation, 'HorizontalAlignment', 'right', 'Fontsize', 11);

linkaxes([ax1 ax2 ax3], 'y')

set(findall(gcf,'-property','FontSize'),'FontSize',12)

axs = findobj(gcf,'type','axe');

width = 0.18;
height = 0.36;
hgap = 0.05;
vgap = 0.03;
voffset = -0.31;

set(axs(3).YLabel,  'Units', 'Normalized', 'Position', [-0.6, 0.5, 0]);

for a = 1:length(axs)
    
    axs(a).Position(3) = width; % width of all subplots
    axs(a).Position(4) = height; % height of all subplots
    
end

% horizontal arrangement
axs(2).Position(1) = axs(3).Position(1) + width + hgap;
axs(1).Position(1) = axs(2).Position(1) + width + hgap;

% vertical arrangement
vpos = axs(3).Position(2) - voffset;
axs(3).Position(2) = vpos;
axs(2).Position(2) = vpos;
axs(1).Position(2) = vpos;

% Align the text labels
hpos = [0.135 0.47 0.81] + 0.03;
for i = 1:3
    set(t1(i), 'Units', 'Normalized')
    t1(i).Position(1:2) = [hpos(i) -0.1];
    set(t2(i), 'Units', 'Normalized')
    t2(i).Position(1:2) = [hpos(i) -0.1];
    set(t3(i), 'Units', 'Normalized')
    t3(i).Position(1:2) = [hpos(i) -0.1];
end

title1.Position(1) = title1.Position(1)-0.2;
title2.Position(1) = title2.Position(1)-0.2;
title3.Position(1) = title3.Position(1)-0.2;

switch sp_dist
    case 'aggregated'
        axs(1).YLim = [-1.2 0.5];
        y1.Position(1) = y1.Position(1)*0.75;
        title1.Position(2) = axs(1).Position(2)+0.16;
        title2.Position(2) = axs(1).Position(2)+0.16;
        title3.Position(2) = axs(1).Position(2)+0.16;
    case 'random'
        axs(1).YLim = [-0.7 0.7];
        y1.Position(1) = y1.Position(1)*0.6;
end

%% Plot Figure S3 (Appendix S2)

plot_boxplots_overview(out_none, out_Ninc, out_Sinc, out_einc)

%% Get results (Figure 4): Scale-dependence of site-selection biases 

% prepare simulation
nvar          = 15; % number of variation steps (between min and max parameter)
output_folder = cd;
job_name      = 'scale_dependency';
job_ID        = '001';
start_coordinates = 'aggregated'; % options: 'random' or 'aggregated'

% start simulation of grain size variation
variable = 'g';
cluster_get_scale_dependency(output_folder, job_name, job_ID,...
    num2str(r_change_rel), num2str(nsamples), variable, num2str(nvar), start_coordinates);
load(sprintf('%s/workspace_%s_%s_%s', output_folder, job_name, job_ID, variable))

% re-scale grain size g to relative percent area
var_none.var = 100*(var_none.var.^2)/(c*c); 
var_Sdec.var = 100*(var_Sdec.var.^2)/(c*c);
var_Sinc.var = 100*(var_Sinc.var.^2)/(c*c);
var_Sinc_g = var_Sinc; % for Figure 5

%% Plot Figure 4: Scale-dependence

plot_variation_horizontal(var_none, var_Sinc, 'g', 100*g*g/(c*c), ...
    'plot_decrease', var_Sdec, 'yticks', [-0.5 0 0.5], 'ylim', [-0.6 0.6])

%% Get results: Additional parameter variations (Figure 5)

% start variation of regional species pool (Stot)
variable = 'Stot';
cluster_get_scale_dependency(output_folder, job_name, job_ID,...
    num2str(r_change_rel), num2str(nsamples), variable, num2str(nvar), start_coordinates);
load(sprintf('%s/workspace_%s_%s_%s', output_folder, job_name, job_ID, variable))
var_Sinc_Stot = var_Sinc;

% start variation of sampling duration (timelag)
variable = 'timelag';
cluster_get_scale_dependency(output_folder, job_name, job_ID,...
    num2str(r_change_rel), num2str(nsamples), variable, num2str(nvar), start_coordinates);
load(sprintf('%s/workspace_%s_%s_%s', output_folder, job_name, job_ID, variable))
var_Sinc_timelag = var_Sinc;
var_Sinc_timelag.var = var_Sinc_timelag.var+1;

% start variation of sampling duration (tend)
variable = 'tend';
cluster_get_scale_dependency(output_folder, job_name, job_ID,...
    num2str(r_change_rel), num2str(nsamples), variable, num2str(nvar), start_coordinates);
load(sprintf('%s/workspace_%s_%s_%s', output_folder, job_name, job_ID, variable))
var_Sinc_tend = var_Sinc;


%% Plot Figure 5 (Impact of bias depends on scale, dispersal & species pool)
% The impact of the bias is calculated in "helperfunc1"

figure('color', 'white', 'position', [217,473,954,134])
subplot(1,7,1)
[ax1] = helperfunc1(var_Sinc_g);
xlabel('Grain size [%]')
subplot(1,7,2)
[ax2] = helperfunc1(var_Sinc_tend);
xlabel('Sampling duration')

subplot(1,7,3)
[ax3] = helperfunc1(var_Sinc_timelag);
ax4.XTick = 2:4:10;
ax4.XLim = [1 10];
xlabel('Dispersal ability')

subplot(1,7,4)
[ax4] = helperfunc1(var_Sinc_Stot);
xlabel('No. species')

linkaxes([ax1, ax2, ax3, ax4], 'y')
ax1.YLim = [0 0.7];
ax1.YTick = [0 0.25 0.5];
ax1.YLabel.String = sprintf('Effect of site-\nselection bias');

% Add triangles
ax = ax1;
xAnnotation = ax.Position(1) + ((100*g*g/(c*c) - ax.XLim(1))/(ax.XLim(2)-ax.XLim(1))) * ax.Position(3);
yAnnotation = ax.Position(2) + ((0 - ax.YLim(1))/(ax.YLim(2)-ax.YLim(1))) * ax.Position(4);
annotation('textbox', [xAnnotation-0.0078, yAnnotation-0.075, 0, 0], 'String', '▲', 'Fontsize', 5, ...
    'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom')
ax = ax2;
xAnnotation = ax.Position(1) + ((tend - ax.XLim(1))/(ax.XLim(2)-ax.XLim(1))) * ax.Position(3);
yAnnotation = ax.Position(2) + ((0 - ax.YLim(1))/(ax.YLim(2)-ax.YLim(1))) * ax.Position(4);
annotation('textbox', [xAnnotation-0.0078, yAnnotation-0.075, 0, 0], 'String', '▲', 'Fontsize', 5, ...
    'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom')
ax = ax3;
xAnnotation = ax.Position(1) + ((1 - ax.XLim(1))/(ax.XLim(2)-ax.XLim(1))) * ax.Position(3);
yAnnotation = ax.Position(2) + ((0 - ax.YLim(1))/(ax.YLim(2)-ax.YLim(1))) * ax.Position(4);
annotation('textbox', [xAnnotation-0.0078, yAnnotation-0.075, 0, 0], 'String', '▲', 'Fontsize', 5, ...
    'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom')
ax = ax4;
xAnnotation = ax.Position(1) + ((Stot - ax.XLim(1))/(ax.XLim(2)-ax.XLim(1))) * ax.Position(3);
yAnnotation = ax.Position(2) + ((0 - ax.YLim(1))/(ax.YLim(2)-ax.YLim(1))) * ax.Position(4);
annotation('textbox', [xAnnotation-0.0078, yAnnotation-0.075, 0, 0], 'String', '▲', 'Fontsize', 5, ...
    'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom')

annotation('textbox', [0.114, 0.999, 0, 0], 'string', 'a', 'Fontsize', 11, 'Fontweight', 'bold', 'Tag', 'subpanel')
annotation('textbox', [0.227, 0.999, 0, 0], 'string', 'b', 'Fontsize', 11, 'Fontweight', 'bold', 'Tag', 'subpanel')
annotation('textbox', [0.344, 0.999, 0, 0], 'string', 'c', 'Fontsize', 11, 'Fontweight', 'bold', 'Tag', 'subpanel')
annotation('textbox', [0.458, 0.999, 0, 0], 'string', 'd', 'Fontsize', 11, 'Fontweight', 'bold', 'Tag', 'subpanel')

l = legend(ax4, 'Rich. bias', 'Comm. bias');
legend('boxoff');
l.Position = [0.5786,0.4291,0.099,0.227];


%% Figure S1
% Experimenting with SAD parameters Stot, mu and var

np = 20; % number of variations
nr = 10000; % number of random repetitions

% Stot
pvec_Stot = round(linspace(50,150,np));
result_Stot = NaN(nr,np);
even_Stot = NaN(nr,np);
SPIE_Stot = NaN(nr,np);
minSAD_Stot = NaN(nr,np);
for i = 1:np
    for r = 1:nr 
       SAD = get_SAD(pvec_Stot(i), 'lognormal', 'mu', mu, 'var', sd);
       result_Stot(r,i) = sum(SAD);
       pi = SAD/sum(SAD);
       PIE   = 1-sum(pi.^2);
       SPIE = 1/(1-PIE); 
       even_Stot(r,i) = PIE;
       SPIE_Stot(r,i) = SPIE;
       minSAD_Stot(r, i) = min(SAD);
    end
end

% mu
pvec_mu = linspace(3,5,np);
result_mu = NaN(nr,np);
even_mu = NaN(nr,np);
SPIE_mu = NaN(nr,np);
minSAD_mu = NaN(nr,np);
for i = 1:np
    for r = 1:nr 
       SAD = get_SAD(Stot, 'lognormal', 'mu', pvec_mu(i), 'var', sd);
       result_mu(r,i) = sum(SAD);
       pi = SAD/sum(SAD);
       PIE   = 1-sum(pi.^2);
       SPIE = 1/(1-PIE); 
       even_mu(r,i) = PIE;
       SPIE_mu(r,i) = SPIE;
       minSAD_mu(r, i) = min(SAD);
    end
end

% var
pvec_var = linspace(1,2,np);
result_var = NaN(nr,np);
even_var = NaN(nr,np);
SPIE_var = NaN(nr,np);
minSAD_var = NaN(nr,np);
for i = 1:np
    for r = 1:nr
       SAD = get_SAD(Stot, 'lognormal', 'mu', mu, 'var', pvec_var(i));
       result_var(r,i) = sum(SAD);
       pi = SAD/sum(SAD);
       PIE   = 1-sum(pi.^2);
       SPIE = 1/(1-PIE); 
       even_var(r,i) = PIE;
       SPIE_var(r,i) = SPIE;
       minSAD_var(r, i) = min(SAD);
    end
end

figure('color', 'white', 'position', [440,205,516,344])
subplot(3,3,1)
ax1 = gca;
plot(pvec_Stot, mean(result_Stot), '.-')
ax1.XTickLabel = [];
axis tight
ylabel(sprintf('Total number\nof individuals'))
subplot(3,3,2)
ax2 = gca;
plot(pvec_mu, mean(result_mu), '.-')
ax2.YTickLabel = []; ax2.XTickLabel = [];
axis tight
subplot(3,3,3)
ax3 = gca;
plot(pvec_var, mean(result_var), '.-')
ax3.YTickLabel = []; ax3.XTickLabel = [];
axis tight
linkaxes([ax1, ax2, ax3], 'y')

subplot(3,3,4)
ax1 = gca;
plot(pvec_Stot, mean(minSAD_Stot), '.-')
ax1.XTickLabel = [];
axis tight
ylabel(sprintf('Abundance of\nrarest species'))
subplot(3,3,5)
ax2 = gca;
plot(pvec_mu, mean(minSAD_mu), '.-')
ax2.YTickLabel = []; ax2.XTickLabel = [];
axis tight
subplot(3,3,6)
ax3 = gca;
plot(pvec_var, mean(minSAD_var), '.-')
ax3.YTickLabel = []; ax3.XTickLabel = [];
axis tight
linkaxes([ax1, ax2, ax3], 'y')

subplot(3,3,7)
ax1 = gca;
plot(pvec_Stot, mean(SPIE_Stot), '.-')
axis tight
xlabel('Total number of species')
ylabel(sprintf('Effective number\nof species (SPIE)'))
subplot(3,3,8)
ax2 = gca;
plot(pvec_mu, mean(SPIE_mu), '.-')
ax2.YTickLabel = [];
axis tight
xlabel('Mean of log abundances')
subplot(3,3,9)
ax3 = gca;
plot(pvec_var, mean(SPIE_var), '.-')
ax3.YTickLabel = [];
axis tight
xlabel('Std of log abundances')
linkaxes([ax1, ax2, ax3], 'y')
annotation('textbox', [0.10, 0.99, 0, 0], 'string', 'a', 'Fontsize', 11, 'Fontweight', 'bold', 'Tag', 'subpanel')
annotation('textbox', [0.38, 0.99, 0, 0], 'string', 'b', 'Fontsize', 11, 'Fontweight', 'bold', 'Tag', 'subpanel')
annotation('textbox', [0.66, 0.99, 0, 0], 'string', 'c', 'Fontsize', 11, 'Fontweight', 'bold', 'Tag', 'subpanel')
annotation('textbox', [0.10, 0.67, 0, 0], 'string', 'd', 'Fontsize', 11, 'Fontweight', 'bold', 'Tag', 'subpanel')
annotation('textbox', [0.38, 0.67, 0, 0], 'string', 'e', 'Fontsize', 11, 'Fontweight', 'bold', 'Tag', 'subpanel')
annotation('textbox', [0.66, 0.67, 0, 0], 'string', 'f', 'Fontsize', 11, 'Fontweight', 'bold', 'Tag', 'subpanel')
annotation('textbox', [0.10, 0.38, 0, 0], 'string', 'g', 'Fontsize', 11, 'Fontweight', 'bold', 'Tag', 'subpanel')
annotation('textbox', [0.38, 0.38, 0, 0], 'string', 'h', 'Fontsize', 11, 'Fontweight', 'bold', 'Tag', 'subpanel')
annotation('textbox', [0.66, 0.38, 0, 0], 'string', 'i', 'Fontsize', 11, 'Fontweight', 'bold', 'Tag', 'subpanel')

%% Figure S2
% Plot number of individuals across grain sizes

gsquare_vec = linspace(0.1, 2, 10);
g_vec = round(sqrt(gsquare_vec*c*c/100));

change_type = 'none';
change_rate_rel = 0;
doAR = 0;
                                   
var_none_g = get_parameter_variation(Stot, c, g, tend, nsamples,...
    change_type, change_rate_rel, g_vec, 'g', doAR, 'mu', mu, 'var', var, 'aggregated');
g_vec_sc = 100*(g_vec.^2)/(c^2)

data = var_none_g;

figure('color', 'white', 'position', [-1461,288,691,147])
subplot(1,3,1)
y = mean(data.N_mean_random');
[lg, pg] = fplot_bounded(data.var.*data.var*100/(c*c), y,...
    y-sd(data.N_mean_random'), y+sd(data.N_mean_random'),...
    mycolors('blue'), 1.5, '-', 0.5);
hold on
y = mean(data.tN_mean_random');
[lt, pt] = fplot_bounded(data.var.*data.var*100/(c*c), y,...
    y-sd(data.tN_mean_random'), y+sd(data.tN_mean_random'),...
    mycolors('green'), 1.5, '-', 0.5);
ax = gca; 
scatter(100*g*g/(c*c), 1, 20, 'k^', 'filled')
set(gca, 'YScale', 'log', 'XScale', 'lin', 'XLim', [-inf inf], ...
    'YLim', [1 inf], 'YTick', [0 1 10 100 1000 10000])
xlabel('Grain size [%]'), ylabel('number of individuals N')
l = legend([lt lg], 'Total', 'Observed in grain', 'Location', 'EastOutside');
l.Position = [0.366275,0.576530,0.1635311,0.18027];


%% Figure S4 and S7: Visualize landscapes aggregated versus Random, versus heterogeneous

% start_coordinates = 'random'; % Figure S7a
% start_coordinates = 'aggregated'; % Figure S7b
start_coordinates = 'heterogeneous'; % Figure S4

nsamples = 1;
change = 'none';
change_rate = 0;
SAD = get_SAD(Stot, 'lognormal', 'mu', mu, 'var', sd, 'visualize', 'off');
[posX, posY, SAD] = simulate_timeseries(SAD, c, g, tend, change, change_rate, start_coordinates);
species_in = repelem(1:Stot, SAD)';
species = repmat(species_in, 1, tend);
visual = get_visual_specifications(species);

figure('color', 'white', 'position', [-654,484,344,305])
hold on
t = 1;  
set(gca, 'XLim', [0 c], 'YLim', [0 c], 'XColor', 'none', 'YColor', 'none')
xlabel(' ')
ylabel(' ')
splot = scatter(posX(:,t), posY(:,t), visual.szs(:,t), visual.colorder(species(:,t),:), 'filled');
splot.MarkerFaceAlpha = 0.8; 
set(gca, 'XTick', [0:g:c], 'YTick', [0:g:c])
grid on
set(gca, 'XTickLabel', [], 'YTickLabel', [], 'GridColor', 'k',...
    'GridAlpha', 1, 'Layer', 'top')
drawnow

%% Figure S9: Heatmap of 4-part-forest

% generate data for Figure S9 a anc c by setting start_coordinates in cell
% above to "random"
% generate data for Figure S9 b anc d by setting start_coordinates in cell
% above to "aggregated"

num_grains = (floor(c/g)^2);
gsteps = (0:((c/g)-1))*g;
[gx_all_mesh, gy_all_mesh] = meshgrid(gsteps, gsteps);
gx_all = gx_all_mesh(:);
gy_all = gy_all_mesh(:);
    
% get N per potential site
comm_t0_N = zeros(num_grains,1);
comm_t0_x = posX(:,1);
comm_t0_y = posY(:,1);
for i = 1:length(comm_t0_x)
    ind = (gx_all == floor(comm_t0_x(i)/g)*g) & ...
    (gy_all == floor(comm_t0_y(i)/g)*g);
    comm_t0_N(ind) = comm_t0_N(ind) + 1;
end

% get S per potential site
visual = get_visual_specifications(species);
rich_t0 = zeros(num_grains,1); % richness per grain
for i = 1:Stot
    occs = unique([floor(posX(species(:,1)==i,1)/g)*g ...
        floor(posY(species(:,1)==i,1)/g)*g], 'rows');
    ind = ismember([gx_all gy_all], occs, 'rows');
    rich_t0(ind) = rich_t0(ind) + 1;
end


figure('color', 'white', 'position', [-820,514,255,319])
subplot(2,1,1)
% histogram(comm_t0_N)
% xlabel('Total abundance at potential site')
% ylabel('# potential sampling sites')
imagesc(reshape(comm_t0_N, c/g, c/g))
set(gca,'YDir','normal');
ax = gca;
ax.Position(3) = 0.7;
colorbar
caxis([9 99])
set(gca,'XTick', [], 'YTick', [])
cb = colorbar;
cb.Label.String = sprintf('Abundance\nper potential site');
subplot(2,1,2)
% histogram(rich_t0)
% xlabel('Richness at potential site')
% ylabel('# potential sampling sites')
imagesc(reshape(rich_t0, c/g, c/g))
set(gca,'YDir','normal');
ax = gca;
ax.Position(3) = 0.7;
cb = colorbar;
cb.Label.String = sprintf('Richness\nper potential site');
caxis([3 44])
set(gca,'XTick', [], 'YTick', [])
set(findall(gcf,'-property','FontSize'),'FontSize',12)

%% Get data for Figure S6

start_coordinates = 'heterogeneous';

% get the effect of bias in this forest
output_folder = cd;
job_name = 'heterogeneous_landscape';
job_ID = 'ID01';
cluster_get_effect_of_biases(output_folder, job_name, job_ID, ...
    num2str(r_change_rel), num2str(nsamples), num2str(c), num2str(g),...
    num2str(Stot), num2str(mu), num2str(sd), num2str(tend),...
    num2str(timelag), 'heterogeneous');
load(sprintf('%s/workspace_%s_%s.mat', output_folder, job_name, job_ID))

%% Plot Figure S6 

fc = [0.5 0.5 0.5;
    mycolors('lightorange');
    mycolors('redorange')];

yl_col = [.8 .8 .8];
yl_width = 2;

boxwidth = 0.1;

trotation = 60;

figure('color', 'white', 'position', [-781,263,459,265])
subplot(1,3,1)
ax1 = gca;
subdata = [out_none.random.S_slope out_none.comm_biased.S_slope out_none.rich_biased.S_slope];
fplot_boxplot(subdata, ax1) 
set(gca, 'XAxisLocation', 'origin', 'Box', 'off')
y1 = ylabel(sprintf('Estimated\nrichness change'));
title1 = title('{\bf a}    No change', 'FontWeight', 'bold', 'HorizontalAlignment', 'left');
t1 = text(ax1, [0.1 0.2 0.3], repmat(ax1.YLim(1)*1.1, 1, 3),...
    {'Random', 'Comm. bias', 'Rich. bias'},...
    'Rotation', trotation, 'HorizontalAlignment', 'right', 'Fontsize', 11);

subplot(1,3,2)
ax2 = gca;
subdata = [out_Sdec.random.S_slope out_Sdec.comm_biased.S_slope out_Sdec.rich_biased.S_slope];
fplot_boxplot(subdata, ax2) 
set(gca, 'XAxisLocation', 'origin', 'Box', 'off', 'XTickLabels', [], 'YTickLabels', [])
title2 = title('{\bf b}   S decrease', 'FontWeight', 'bold', 'HorizontalAlignment', 'left');
t2 = text(ax2, [0.1 0.2 0.3], repmat(ax2.YLim(1)*1.1, 1, 3),...
    {'Random', 'Comm. bias', 'Rich. bias'},...
    'Rotation', trotation, 'HorizontalAlignment', 'right', 'Fontsize', 11);

subplot(1,3,3)
ax3 = gca;
subdata = [out_Sinc.random.S_slope out_Sinc.comm_biased.S_slope out_Sinc.rich_biased.S_slope];
fplot_boxplot(subdata, ax3) 
set(gca, 'XAxisLocation', 'origin', 'Box', 'off', 'XTickLabels', [], 'YTickLabels', [])
title3 = title('{\bf c}    S increase', 'FontWeight', 'bold', 'HorizontalAlignment', 'left');
t3 = text(ax3, [0.1 0.2 0.3], repmat(ax3.YLim(1)*1.1, 1, 3),...
    {'Random', 'Comm. bias', 'Rich. bias'},...
    'Rotation', trotation, 'HorizontalAlignment', 'right', 'Fontsize', 11);

linkaxes([ax1 ax2 ax3], 'y')

set(findall(gcf,'-property','FontSize'),'FontSize',12)

axs = findobj(gcf,'type','axe');

width = 0.18;
height = 0.36;
hgap = 0.05;
vgap = 0.03; %0.19;
voffset = -0.31;

set(axs(3).YLabel,  'Units', 'Normalized', 'Position', [-0.6, 0.5, 0]);

for a = 1:length(axs)
    
    axs(a).Position(3) = width; % width of all subplots
    axs(a).Position(4) = height; % height of all subplots
    
end

% horizontal arrangement
axs(2).Position(1) = axs(3).Position(1) + width + hgap;
axs(1).Position(1) = axs(2).Position(1) + width + hgap;

% vertical arrangement
vpos = axs(3).Position(2) - voffset;
axs(3).Position(2) = vpos;
axs(2).Position(2) = vpos;
axs(1).Position(2) = vpos;

% Align the text labels
hpos = [0.135 0.47 0.81] + 0.03;
for i = 1:3
    set(t1(i), 'Units', 'Normalized')
    t1(i).Position(1:2) = [hpos(i) -0.1];
    set(t2(i), 'Units', 'Normalized')
    t2(i).Position(1:2) = [hpos(i) -0.1];
    set(t3(i), 'Units', 'Normalized')
    t3(i).Position(1:2) = [hpos(i) -0.1];
end

title1.Position(1) = title1.Position(1)-0.2;
title2.Position(1) = title2.Position(1)-0.2;
title3.Position(1) = title3.Position(1)-0.2;

axs(1).YLim = [-0.7 0.7];
y1.Position(1) = y1.Position(1)*0.6;


%% Figure S5
% Get effect of bias WITHIN each quadrant ("specialist" case)

% select a quadrant
quadrant = 4;

c  = 50; % the quadrant is 50x50 (default landscape = 100x100)

if quadrant == 1
    mu = 3.5;
    sd = 1.5;
    Stot = 100;
elseif quadrant == 2
    mu = 1.2;
    sd = 1.5;
    Stot = 100;
elseif quadrant == 3
    mu = 3.5;
    sd = 1.5;
    Stot = 50;
elseif quadrant == 4
    mu = 1.2;
    sd = 1.5;
    Stot = 50;
end

% Start simulation (do this for each quadrant)
output_folder = cd;
job_name      = 'effect_of_biases';
job_ID        = '001';
cluster_get_effect_of_biases(output_folder, job_name, job_ID, ...
    num2str(r_change_rel), num2str(nsamples), num2str(c), num2str(g),...
    num2str(Stot), num2str(mu), num2str(sd), num2str(tend),...
    num2str(timelag), start_coordinates);
load(sprintf('%s/workspace_%s_%s.mat', output_folder, job_name, job_ID))
   
% Plot (one row of Figure S5):       
fc = [0.5 0.5 0.5;
    mycolors('lightorange');
    mycolors('redorange')];

yl_col = [.8 .8 .8];
yl_width = 2;

boxwidth = 0.1;

trotation = 60;

figure('color', 'white', 'position', [-781,263,459,265])
subplot(1,3,1)
ax1 = gca;
subdata = [out_none.random.S_slope out_none.comm_biased.S_slope out_none.rich_biased.S_slope];
fplot_boxplot(subdata, ax1) 
set(gca, 'XAxisLocation', 'origin', 'Box', 'off')
y1 = ylabel(sprintf('Estimated\nrichness change'));

subplot(1,3,2)
ax2 = gca;
subdata = [out_Sdec.random.S_slope out_Sdec.comm_biased.S_slope out_Sdec.rich_biased.S_slope];
fplot_boxplot(subdata, ax2) 
set(gca, 'XAxisLocation', 'origin', 'Box', 'off', 'XTickLabels', [], 'YTickLabels', [])

subplot(1,3,3)
ax3 = gca;
subdata = [out_Sinc.random.S_slope out_Sinc.comm_biased.S_slope out_Sinc.rich_biased.S_slope];
fplot_boxplot(subdata, ax3) 
set(gca, 'XAxisLocation', 'origin', 'Box', 'off', 'XTickLabels', [], 'YTickLabels', [])

linkaxes([ax1 ax2 ax3], 'y')

set(findall(gcf,'-property','FontSize'),'FontSize',12)

axs = findobj(gcf,'type','axe');

width = 0.18;
height = 0.36;
hgap = 0.05;
vgap = 0.03; %0.19;
voffset = -0.31;

set(axs(3).YLabel,  'Units', 'Normalized', 'Position', [-0.6, 0.5, 0]);

for a = 1:length(axs)
    
    axs(a).Position(3) = width; % width of all subplots
    axs(a).Position(4) = height; % height of all subplots
    
end

% horizontal arrangement
axs(2).Position(1) = axs(3).Position(1) + width + hgap;
axs(1).Position(1) = axs(2).Position(1) + width + hgap;

% vertical arrangement
vpos = axs(3).Position(2) - voffset;
axs(3).Position(2) = vpos;
axs(2).Position(2) = vpos;
axs(1).Position(2) = vpos;

axs(1).YLim = [-0.7 0.7];
y1.Position(1) = y1.Position(1)*0.6;








