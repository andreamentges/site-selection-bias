%% Plot overview of boxplots 
% Plots boxplots for all diversity change regimes (none, richness, and
% evenness), all diversity measures (N, S, SPIE), and for the three
% sampling strategies (random, community-, and richness-biased)

function plot_boxplots_overview(out_none, out_N, out_S, out_evenness)
%% Boxplots

fc = [0.5 0.5 0.5;
    mycolors('lightorange');
    mycolors('redorange')];

yl_col = [.8 .8 .8];
yl_width = 2;

boxwidth = 0.1;

%%
%%%%%%%%%%%
% S
%%%%%%%%%%%
figure('color', 'white', 'position', [440,78,459,778])
subplot(4,4,1)
ax1 = gca;
subdata = [out_none.random.S_slope out_none.comm_biased.S_slope out_none.rich_biased.S_slope];
plot([ax1.XLim(1) ax1.XLim(2)], [median(subdata(:,1)) median(subdata(:,1))], 'color', yl_col, 'linewidth', yl_width)
hold on
fplot_boxplot(subdata, ax1) 
set(gca, 'XAxisLocation', 'origin', 'Box', 'off')
ylabel('slope S')
title('No change')

subplot(4,4,2)
ax2 = gca;
subdata = [out_S.random.S_slope out_S.comm_biased.S_slope out_S.rich_biased.S_slope];
plot([ax1.XLim(1) ax1.XLim(2)], [median(subdata(:,1)) median(subdata(:,1))], 'color', yl_col, 'linewidth', yl_width)
hold on
fplot_boxplot(subdata, ax2) 
set(gca, 'XAxisLocation', 'origin', 'Box', 'off', 'XTickLabels', [], 'YTickLabels', [])
title('S increase')


subplot(4,4,3)
ax3 = gca;
subdata = [out_N.random.S_slope out_N.comm_biased.S_slope out_N.rich_biased.S_slope];
plot([ax1.XLim(1) ax1.XLim(2)], [median(subdata(:,1)) median(subdata(:,1))], 'color', yl_col, 'linewidth', yl_width)
hold on
fplot_boxplot(subdata, ax3) 
set(gca, 'XAxisLocation', 'origin', 'Box', 'off', 'XTickLabels', [], 'YTickLabels', [])
title('N increase')

subplot(4,4,4)
ax4 = gca;
subdata = [out_evenness.random.S_slope out_evenness.comm_biased.S_slope out_evenness.rich_biased.S_slope];
plot([ax1.XLim(1) ax1.XLim(2)], [median(subdata(:,1)) median(subdata(:,1))], 'color', yl_col, 'linewidth', yl_width)
hold on
fplot_boxplot(subdata, ax4) 
set(gca, 'XAxisLocation', 'origin', 'Box', 'off', 'XTickLabels', [], 'YTickLabels', [])
title('Evenn. increase')

linkaxes([ax1 ax2 ax3 ax4], 'y'), set(gca, 'YLim', [-1.1 1])
ax1.YTick = [-1 -0.5 0 0.5];

%%%%%%%%
% N
%%%%%%%%
subplot(4,4,5)
ax5 = gca;
subdata = [out_none.random.N_slope out_none.comm_biased.N_slope out_none.rich_biased.N_slope];
plot([ax1.XLim(1) ax1.XLim(2)], [median(subdata(:,1)) median(subdata(:,1))], 'color', yl_col, 'linewidth', yl_width)
hold on
fplot_boxplot(subdata, ax5)  
set(gca, 'XAxisLocation', 'origin', 'Box', 'off', 'XTickLabels', [])
ylabel('slope N')

subplot(4,4,6)
ax6 = gca;
subdata = [out_S.random.N_slope out_S.comm_biased.N_slope out_S.rich_biased.N_slope];
plot([ax1.XLim(1) ax1.XLim(2)], [median(subdata(:,1)) median(subdata(:,1))], 'color', yl_col, 'linewidth', yl_width)
hold on
fplot_boxplot(subdata, ax6) 
set(gca, 'XAxisLocation', 'origin', 'Box', 'off', 'XTickLabels', [], 'YTickLabels', [])

subplot(4,4,7)
ax7 = gca;
subdata = [out_N.random.N_slope out_N.comm_biased.N_slope out_N.rich_biased.N_slope];
plot([ax1.XLim(1) ax1.XLim(2)], [median(subdata(:,1)) median(subdata(:,1))], 'color', yl_col, 'linewidth', yl_width)
hold on
fplot_boxplot(subdata, ax7) 
set(gca, 'XAxisLocation', 'origin', 'Box', 'off', 'XTickLabels', [], 'YTickLabels', [])

subplot(4,4,8)
ax8 = gca;
subdata = [out_evenness.random.N_slope out_evenness.comm_biased.N_slope out_evenness.rich_biased.N_slope];
plot([ax1.XLim(1) ax1.XLim(2)], [median(subdata(:,1)) median(subdata(:,1))], 'color', yl_col, 'linewidth', yl_width)
hold on
fplot_boxplot(subdata, ax8) 
set(gca, 'XAxisLocation', 'origin', 'Box', 'off', 'XTickLabels', [], 'YTickLabels', [])

linkaxes([ax5 ax6 ax7 ax8], 'y'), set(gca, 'YLim', [-2.2 1.5])

%%%%%%%%
% SPIE
%%%%%%%%
subplot(4,4,9)
ax9 = gca;
subdata = [out_none.random.S_PIE_slope out_none.comm_biased.S_PIE_slope out_none.rich_biased.S_PIE_slope];
plot([ax1.XLim(1) ax1.XLim(2)], [median(subdata(:,1)) median(subdata(:,1))], 'color', yl_col, 'linewidth', yl_width)
hold on
fplot_boxplot(subdata, ax9);% 'labels', {'Random', 'Comm. bias', 'Rich. bias'})
set(gca, 'XAxisLocation', 'origin', 'Box', 'off', 'XTickLabels', [])
ylabel('slope SPIE', 'Interpreter', 'none')
t1 = text(ax9, [0.1 0.2 0.3], repmat(ax9.YLim(1)*1.1, 1, 3),...
    {'Random', 'Comm. bias', 'Rich. bias'},...
    'Rotation', 90, 'HorizontalAlignment', 'right', 'Fontsize', 11);

subplot(4,4,10)
ax10 = gca;
subdata = [out_S.random.S_PIE_slope out_S.comm_biased.S_PIE_slope out_S.rich_biased.S_PIE_slope];
plot([ax1.XLim(1) ax1.XLim(2)], [median(subdata(:,1)) median(subdata(:,1))], 'color', yl_col, 'linewidth', yl_width)
hold on
fplot_boxplot(subdata, ax10)
set(gca, 'XAxisLocation', 'origin', 'Box', 'off', 'XTickLabels', [], 'YTickLabels', [])
t2 = text(ax10, [0.1 0.2 0.3], repmat(ax10.YLim(1)*1.1, 1, 3),...
    {'Random', 'Comm. bias', 'Rich. bias'},...
    'Rotation', 90, 'HorizontalAlignment', 'right', 'Fontsize', 11);

subplot(4,4,11)
ax11 = gca;
subdata = [out_N.random.S_PIE_slope out_N.comm_biased.S_PIE_slope out_N.rich_biased.S_PIE_slope];
plot([ax1.XLim(1) ax1.XLim(2)], [median(subdata(:,1)) median(subdata(:,1))], 'color', yl_col, 'linewidth', yl_width)
hold on
fplot_boxplot(subdata, ax11)
set(gca, 'XAxisLocation', 'origin', 'Box', 'off', 'XTickLabels', [], 'YTickLabels', [])
t3 = text(ax11, [0.1 0.2 0.3], repmat(ax11.YLim(1)*1.1, 1, 3),...
    {'Random', 'Comm. bias', 'Rich. bias'},...
    'Rotation', 90, 'HorizontalAlignment', 'right', 'Fontsize', 11);

subplot(4,4,12)
ax12 = gca;
subdata = [out_evenness.random.S_PIE_slope out_evenness.comm_biased.S_PIE_slope out_evenness.rich_biased.S_PIE_slope];
plot([ax1.XLim(1) ax1.XLim(2)], [median(subdata(:,1)) median(subdata(:,1))], 'color', yl_col, 'linewidth', yl_width)
hold on
fplot_boxplot(subdata, ax12)
set(gca, 'XAxisLocation', 'origin', 'Box', 'off', 'XTickLabels', [], 'YTickLabels', [])
t4 = text(ax12, [0.1 0.2 0.3], repmat(ax12.YLim(1)*1.1, 1, 3),...
    {'Random', 'Comm. bias', 'Rich. bias'},...
    'Rotation', 90, 'HorizontalAlignment', 'right', 'Fontsize', 11);

linkaxes([ax9 ax10 ax11 ax12], 'y'), set(gca, 'YLim', [-1.1 1])
ax9.YTick = [-1 -0.5 0 0.5];


%% Subplot format 

set(findall(gcf,'-property','FontSize'),'FontSize',12)

axs = findobj(gcf,'type','axe');

width = 0.15;
height = 0.17;
gap = 0.07;
vgap = 0.05; %0.19;
voffset = 0.11;

set(axs(12).YLabel, 'Units', 'Normalized', 'Position', [-0.6, 0.5, 0]);
set(axs(8).YLabel,  'Units', 'Normalized', 'Position', [-0.6, 0.5, 0]);
set(axs(4).YLabel,  'Units', 'Normalized', 'Position', [-0.6, 0.5, 0]);

for a = 1:length(axs)
    
    axs(a).Position(3) = width; % width of all subplots
    axs(a).Position(4) = height; % height of all subplots
    
end

% horizontal arrangement
l_c1 = axs(4).Position(1);
axs(12).Position(1) = l_c1;
axs(8).Position(1)  = l_c1;

l_c2 = l_c1 + width + gap;
axs(11).Position(1) = l_c2;
axs(7).Position(1)  = l_c2;
axs(3).Position(1)  = l_c2;

l_c3 = l_c2 + width + gap;
axs(10).Position(1) = l_c3;
axs(6).Position(1)  = l_c3;
axs(2).Position(1)  = l_c3;

l_c4 = l_c3 + width + gap;
axs(9).Position(1) = l_c4;
axs(5).Position(1)  = l_c4;
axs(1).Position(1)  = l_c4;

% vertical arrangement
d_r1 = axs(12).Position(2) - voffset;
axs(12).Position(2) = d_r1;
axs(11).Position(2) = d_r1;
axs(10).Position(2) = d_r1;
axs(9).Position(2)  = d_r1;

d_r2 = axs(12).Position(2) - height - vgap;
axs(8).Position(2) = d_r2;
axs(7).Position(2) = d_r2;
axs(6).Position(2) = d_r2;
axs(5).Position(2) = d_r2;

d_r3 = d_r2 - height - vgap;
axs(4).Position(2) = d_r3;
axs(3).Position(2) = d_r3;
axs(2).Position(2) = d_r3;
axs(1).Position(2) = d_r3;

% Align the text labels
hpos = [0.135 0.47 0.81];
for i = 1:3
    set(t1(i), 'Units', 'Normalized')
    t1(i).Position(1:2) = [hpos(i) -0.1];
    set(t2(i), 'Units', 'Normalized')
    t2(i).Position(1:2) = [hpos(i) -0.1];
    set(t3(i), 'Units', 'Normalized')
    t3(i).Position(1:2) = [hpos(i) -0.1];
    set(t4(i), 'Units', 'Normalized')
    t4(i).Position(1:2) = [hpos(i) -0.1];
end

offs = 0.03;
annotation('textbox', [0.1,  0.884-offs, 0, 0], 'string', 'a', 'Fontsize', 11, 'Fontweight', 'bold', 'Tag', 'subpanel')
annotation('textbox', [0.32, 0.884-offs, 0, 0], 'string', 'b', 'Fontsize', 11, 'Fontweight', 'bold', 'Tag', 'subpanel')
annotation('textbox', [0.54, 0.884-offs, 0, 0], 'string', 'c', 'Fontsize', 11, 'Fontweight', 'bold', 'Tag', 'subpanel')
annotation('textbox', [0.74, 0.884-offs, 0, 0], 'string', 'd', 'Fontsize', 11, 'Fontweight', 'bold', 'Tag', 'subpanel')
annotation('textbox', [0.1,  0.66-offs, 0, 0], 'string', 'e', 'Fontsize', 11, 'Fontweight', 'bold', 'Tag', 'subpanel')
annotation('textbox', [0.32, 0.66-offs, 0, 0], 'string', 'f', 'Fontsize', 11, 'Fontweight', 'bold', 'Tag', 'subpanel')
annotation('textbox', [0.54, 0.66-offs, 0, 0], 'string', 'g', 'Fontsize', 11, 'Fontweight', 'bold', 'Tag', 'subpanel')
annotation('textbox', [0.74, 0.66-offs, 0, 0], 'string', 'h', 'Fontsize', 11, 'Fontweight', 'bold', 'Tag', 'subpanel')
annotation('textbox', [0.1,  0.44-offs, 0, 0], 'string', 'i', 'Fontsize', 11, 'Fontweight', 'bold', 'Tag', 'subpanel')
annotation('textbox', [0.32, 0.44-offs, 0, 0], 'string', 'j', 'Fontsize', 11, 'Fontweight', 'bold', 'Tag', 'subpanel')
annotation('textbox', [0.54, 0.44-offs, 0, 0], 'string', 'k', 'Fontsize', 11, 'Fontweight', 'bold', 'Tag', 'subpanel')
annotation('textbox', [0.74, 0.44-offs, 0, 0], 'string', 'l', 'Fontsize', 11, 'Fontweight', 'bold', 'Tag', 'subpanel')



end
