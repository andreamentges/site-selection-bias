%% Draw half-transparent area around line
%
% Plots a line with half-transparent shaded area around it. 
%
% lh, ph, axh are handles to patch, line and axis

function [lh, ph, axh] = plot_bounded(x, y, errmin, errmax, color, linewidth, linestyle, alpha);

ph = patch([x fliplr(x)],...
    [errmax fliplr(errmin)],...
    color, 'EdgeColor', 'none', 'FaceAlpha', alpha);

[~, uniqueIdx] = unique(x); 
    
hold on
lh = plot(x(uniqueIdx'), y(uniqueIdx'), 'color', color, 'linewidth', linewidth, 'linestyle', linestyle);
lh.Color(4) = 0.8;

axh = gca;


end