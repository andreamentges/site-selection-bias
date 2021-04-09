%% Plot the estimated richness change over a gradient of one parameter 
% eg grain size, dispersal ability ... 

function plot_variation_horizontal(var_none, var_Sinc, x_label, default, varargin)

% Optional input:
%
% 'yticks', [-0.4 0 0.4]

%% Plotting parameters

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

plot_decrease = 'off';
if any(strcmp(varargin, 'plot_decrease'))
    plot_decrease = 'on';
    idx = find(strcmp(varargin, 'plot_decrease'));
    var_Sdec = varargin{idx+1};
end 
    

%% NO CHANGE

hoffset = 0.01;

figure('color', 'white', 'position', [42,448,576,149])
if ~isempty(var_none)
    data = var_none;
    if strcmp(x_label, 'timelag')
        % convert time lag to "max. movement speed (cells/timepoint)
        data.var = data.var+1;
    end
    if strcmp(plot_decrease, 'on')
        subplot(1,3,1)
    else
        subplot(1,2,1)
    end
    subdata  = data.S_slope_rich_biased';
    errmin = mean(subdata) + std(subdata);
    errmax = mean(subdata) - std(subdata);
    fplot_bounded(data.var, mean(subdata),...
        errmin, errmax, linecol_obs_r, linewidth_obs_r, linestyle_obs_r, alphar)
    subdata  = data.S_slope_comm_biased';
    errmin = mean(subdata) + std(subdata);
    errmax = mean(subdata) - std(subdata);
    fplot_bounded(data.var, mean(subdata),...
        errmin, errmax, linecol_obs_c, linewidth_obs_c, linestyle_obs_c, alphac)
    subdata  = data.S_slope_random';
    errmin = mean(subdata) + std(subdata);
    errmax = mean(subdata) - std(subdata);
    fplot_bounded(data.var, mean(subdata),...
        errmin, errmax, linecol_obs_n, linewidth_obs_n, linestyle_obs_n, alphan)
    ax_old_A = gca; 
    ax_old_A.XAxisLocation = 'origin'; 
    xticks = ax_old_A.XTick;
    ax_old_A.XLim = [min(data.var) max(data.var)];
    ax_old_A.XTick = [];
    ax_old_A.Position = ax_old_A.Position + [hoffset 0 -hoffset 0];
    ax_new_A = axes('position', ax_old_A.Position, 'color', 'none');
    ax_new_A.XLim = [min(data.var) max(data.var)];
    ax_new_A.XTick = xticks;
    ax_new_A.YLim = ax_old_A.YLim;
    ax_new_A.YTick = ax_old_A.YTick;
    ax_old_A.YTick = [];
    yl1 = ylabel(sprintf('Estimated\nrichness change'));
    xlabel(x_label, 'interpreter', 'none')
        
    ax_new_A.XRuler.Axle.LineStyle = 'none';
    ax_new_A.XRuler.Axle.LineStyle = 'none';
    ax_new_A.XRuler.TickLength = [0.02 0.02];
    
end
tA = title('No change');

if strcmp(plot_decrease, 'on')
    data = var_Sdec;
    if strcmp(x_label, 'timelag')
        % convert time lag to "max. movement speed (cells/timepoint)
        data.var = data.var+1;
    end
    subplot(1,3,2)
    subdata  = data.S_slope_rich_biased';
    errmin = mean(subdata) + std(subdata);
    errmax = mean(subdata) - std(subdata);
    fplot_bounded(data.var, mean(subdata),...
        errmin, errmax, linecol_obs_r, linewidth_obs_r, linestyle_obs_r, alphar)
    subdata  = data.S_slope_comm_biased';
    errmin = mean(subdata) + std(subdata);
    errmax = mean(subdata) - std(subdata);
    fplot_bounded(data.var, mean(subdata),...
        errmin, errmax, linecol_obs_c, linewidth_obs_c, linestyle_obs_c, alphac)
    subdata  = data.S_slope_random';
    errmin = mean(subdata) + std(subdata);
    errmax = mean(subdata) - std(subdata);
    fplot_bounded(data.var, mean(subdata),...
        errmin, errmax, linecol_obs_n, linewidth_obs_n, linestyle_obs_n, alphan)
    ax_old_C = gca; 
    ax_old_C.XAxisLocation = 'origin'; 
    xticks = ax_old_C.XTick;
    ax_old_C.XLim = [min(data.var) max(data.var)];
    ax_old_C.XTick = [];
    ax_old_C.Position = ax_old_C.Position + [hoffset 0 -hoffset 0];
    ax_new_C = axes('position', ax_old_C.Position, 'color', 'none');
    ax_new_C.XLim = [min(data.var) max(data.var)];
    ax_new_C.XTick = xticks;
    ax_new_C.YLim = ax_old_C.YLim;
    ax_new_C.YTick = [];
    ax_old_C.YTick = [];
    xlabel(x_label, 'interpreter', 'none')
        
    ax_new_C.XRuler.Axle.LineStyle = 'none';
    ax_new_C.XRuler.Axle.LineStyle = 'none';
    ax_new_C.XRuler.TickLength = [0.02 0.02];
    
    tC = title('S decrease');
end


% INCREASE
if ~isempty(var_Sinc)
    
    if strcmp(plot_decrease, 'on')
        subplot(1,3,3)
    else
        subplot(1,2,2)
    end
    data = var_Sinc;
    if strcmp(x_label, 'timelag')
        % convert time lag to "max. movement speed (cells/timepoint)
        data.var = data.var+1;
    end
    subdata  = data.S_slope_rich_biased';
    errmin = mean(subdata) + std(subdata);
    errmax = mean(subdata) - std(subdata);
    [lhr] = fplot_bounded(data.var, mean(subdata),...
        errmin, errmax, linecol_obs_r, linewidth_obs_r, linestyle_obs_r, alphar);
    subdata  = data.S_slope_comm_biased';
    errmin = mean(subdata) + std(subdata);
    errmax = mean(subdata) - std(subdata);
    [lhc] = fplot_bounded(data.var, mean(subdata),...
        errmin, errmax, linecol_obs_c, linewidth_obs_c, linestyle_obs_c, alphac);
    subdata  = data.S_slope_random';
    errmin = mean(subdata) + std(subdata);
    errmax = mean(subdata) - std(subdata);
    [lhn] = fplot_bounded(data.var, mean(subdata),...
        errmin, errmax, linecol_obs_n, linewidth_obs_n, linestyle_obs_n, alphan);
    l = legend([lhn, lhc, lhr], 'Random', 'Comm. bias', 'Rich. bias', 'location', 'eastoutside');
    legend('boxoff');
    l.FontSize = 10;
    l.Position = [0.729668531788526,0.348398789733794,0.166666666666666,0.295302013422819];
    ax_old_B = gca; ax_old_B.XAxisLocation = 'origin'; 
    xticks = ax_old_B.XTick;
    ax_old_B.XLim = [min(data.var) max(data.var)];
    ax_old_B.XTick = [];
    ax_old_B.Position = ax_old_B.Position + [hoffset 0 -hoffset 0];
    ax_new_B = axes('position', ax_old_B.Position, 'color', 'none');
    ax_new_B.XLim = [min(data.var) max(data.var)];
    ax_new_B.XTick = xticks;
    ax_new_B.YLim = ax_old_B.YLim;
    ax_new_B.YTick = ax_old_B.YTick;
    
    ax_old_B.YTick = [];
    xlabel(x_label, 'interpreter', 'none')
    ax_new_B.XRuler.Axle.LineStyle = 'none';
    ax_new_B.XRuler.TickLength = [0.02 0.02];
    
end
tB = title('S increase');

if strcmp(plot_decrease, 'on')
    linkaxes([ax_old_A, ax_new_A, ax_old_C, ax_new_C, ax_old_B, ax_new_B], 'y')
else
    if ~isempty(var_none) && ~isempty(var_Sinc)
        linkaxes([ax_old_A, ax_new_A, ax_old_B, ax_new_B], 'y')
    elseif isempty(var_none) && ~isempty(var_Sinc)
        linkaxes([ax_old_B, ax_new_B], 'y')
    end
end

if any(strcmp(varargin, 'yticks'))
    idx = find(strcmp(varargin, 'yticks'));
    yticks = varargin{idx+1};
    ax_new_A.YTick = yticks;
    ax_new_B.YTick = yticks;
    if strcmp(plot_decrease, 'on')
        ax_new_B.YTick = [];
    end
end

if any(strcmp(varargin, 'ylim'))
    idx = find(strcmp(varargin, 'ylim'));
    ylim = varargin{idx+1};
    ax_new_A.YLim = ylim;
    ax_new_B.YLim = ylim;
    if strcmp(plot_decrease, 'on')
       ax_new_C.YLim = ylim; 
    end
end

if strcmp(x_label, 'timelag')
    ax_new_A.XTick = 2:2:10;
    ax_new_B.XTick = 2:2:10;
    ax_new_A.XLim = [1 10];
    ax_new_B.XLim = [1 10];
    ax_new_A.XLabel.String = 'Dispersal ability';
    ax_new_B.XLabel.String = 'Dispersal ability';
    if strcmp(plot_decrease, 'on')
    	ax_new_C.XTick = 2:2:10;
        ax_new_C.XLim = [1 10];
        ax_new_C.XLabel.String = 'Dispersal ability';
    end
end

if strcmp(x_label, 'g')
    ax_new_A.XLabel.String = 'grain size (% of total area)';
    ax_new_B.XLabel.String = 'grain size (% of total area)';
    if strcmp(plot_decrease, 'on')
        ax_new_C.XLabel.String = 'grain size (% of total area)';
        ax_new_A.XLabel.String = '';
        ax_new_B.XLabel.String = '';
    end
end

if strcmp(x_label, 'c')
    ax_new_A.XLabel.String ='landscape size';
    ax_new_B.XLabel.String ='landscape size';
    if strcmp(plot_decrease, 'on')
        ax_new_C.XLabel.String ='landscape size';
    end
end

drawnow

% add triangle marking the default value
if ~isempty(var_none)
    hold(ax_new_A, 'on')
    if strcmp(x_label, 'timelag')
        scatter(ax_new_A, default+1, ax_new_A.YLim(1), 25, 'k^', 'filled')
    else
        scatter(ax_new_A, default, ax_new_A.YLim(1), 25, 'k^', 'filled')
    end
end
if ~isempty(var_Sinc)
    hold(ax_new_B, 'on')
    if strcmp(x_label, 'timelag')
        scatter(ax_new_B, default+1, ax_new_B.YLim(1), 25, 'k^', 'filled')
    else
        scatter(ax_new_B, default, ax_new_B.YLim(1), 25, 'k^', 'filled')
    end
end
if ~isempty(var_Sinc) && ~isempty(var_Sinc) && strcmp(plot_decrease, 'on')
    hold(ax_new_C, 'on')
    if strcmp(x_label, 'timelag')
        scatter(ax_new_C, default+1, ax_new_C.YLim(1), 25, 'k^', 'filled')
    else
        scatter(ax_new_C, default, ax_new_C.YLim(1), 25, 'k^', 'filled')
    end
end

%%

if strcmp(plot_decrease, 'off')
    height = 0.6;
    width = 0.25;
    ypos = 0.22;
    xpos = 0.1;

    ax_old_A.Position(4) = height;
    ax_new_A.Position(4) = height;
    ax_old_B.Position(4) = height;
    ax_new_B.Position(4) = height;

    ax_old_A.Position(3) = width;
    ax_new_A.Position(3) = width;
    ax_old_B.Position(3) = width;
    ax_new_B.Position(3) = width;

    ax_old_A.Position(2) = ypos;
    ax_new_A.Position(2) = ypos;
    ax_old_B.Position(2) = ypos;
    ax_new_B.Position(2) = ypos;

    ax_old_A.Position(1) = xpos;
    ax_new_A.Position(1) = xpos;
    ax_old_B.Position(1) = xpos+0.35;
    ax_new_B.Position(1) = xpos+0.35;
    
    drawnow

    set(findall(gcf,'-property','FontSize'),'FontSize',11)
    l.FontSize = 10;

    annotation('textbox', [0.08, 0.96, 0, 0], 'string', 'a', 'Fontsize', 11, 'Fontweight', 'bold', 'Tag', 'subpanel')
    annotation('textbox', [0.43, 0.96, 0, 0], 'string', 'b', 'Fontsize', 11, 'Fontweight', 'bold', 'Tag', 'subpanel')


else
    
    height = 0.6;
    width = 0.2;
    ypos = 0.23;
    xpos = 0.1;

    ax_old_A.Position(4) = height;
    ax_new_A.Position(4) = height;
    ax_old_B.Position(4) = height;
    ax_new_B.Position(4) = height;
    ax_old_C.Position(4) = height;
    ax_new_C.Position(4) = height;

    ax_old_A.Position(3) = width;
    ax_new_A.Position(3) = width;
    ax_old_B.Position(3) = width;
    ax_new_B.Position(3) = width;
    ax_old_C.Position(3) = width;
    ax_new_C.Position(3) = width;

    ax_old_A.Position(2) = ypos;
    ax_new_A.Position(2) = ypos;
    ax_old_B.Position(2) = ypos;
    ax_new_B.Position(2) = ypos;
    ax_old_C.Position(2) = ypos;
    ax_new_C.Position(2) = ypos;

    ax_old_A.Position(1) = xpos;
    ax_new_A.Position(1) = xpos;
    ax_old_C.Position(1) = xpos+0.23;
    ax_new_C.Position(1) = xpos+0.23;
    ax_old_B.Position(1) = xpos+0.23*2;
    ax_new_B.Position(1) = xpos+0.23*2;
    
    drawnow

    set(findall(gcf,'-property','FontSize'),'FontSize',11)
    l.FontSize = 10;

    annotation('textbox', [0.08, 0.96, 0, 0], 'string', 'a', 'Fontsize', 11, 'Fontweight', 'bold', 'Tag', 'subpanel')
    annotation('textbox', [0.31, 0.96, 0, 0], 'string', 'b', 'Fontsize', 11, 'Fontweight', 'bold', 'Tag', 'subpanel')
    annotation('textbox', [0.54, 0.96, 0, 0], 'string', 'c', 'Fontsize', 11, 'Fontweight', 'bold', 'Tag', 'subpanel')
    
    l.Position = [0.780015754010748,0.381955836713660,0.166666666666666,0.295302013422819];
    
end





end

