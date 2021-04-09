%% Helper function is for plotting Figure 5 
% (impact of bias depends on scale, dispersal, and species pool size)

function [ax] = helperfunc1(data);


linewidth_obs_n = 2.5;
linewidth_obs_c = 1.5;
linewidth_obs_r = 1;

linecol_obs_n = [.6 .6 .6];
linecol_obs_c = mycolors('lightorange');
linecol_obs_r = mycolors('red');

linestyle_obs_n =  ':';
linestyle_obs_c =  '-';
linestyle_obs_r =  '-';

alphar = 0.10;
alphac = 0.15;
alphan = 0.15;

subdata_n  = data.S_slope_random;

subdata  =  subdata_n' - data.S_slope_rich_biased';
lh = plot(data.var, mean(subdata), 'color', linecol_obs_r, 'linewidth', linewidth_obs_r, 'linestyle', linestyle_obs_r);
lh.Color(4) = 0.8;
hold on

subdata  = subdata_n' - data.S_slope_comm_biased';
lh = plot(data.var, mean(subdata), 'color', linecol_obs_c, 'linewidth', linewidth_obs_c, 'linestyle', linestyle_obs_c);
lh.Color(4) = 0.8;
    

ax = gca;
    
ax.XLim = [min(data.var) max(data.var)];
ax.YTick = [];
ax.XAxisLocation = 'bottom';
ax.Position(2) = ax.Position(2)+0.14;
ax.Position(4) = ax.Position(4)-0.2;


    
end