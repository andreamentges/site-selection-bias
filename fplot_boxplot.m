%% Draw a boxplot with additional points and horizontal bar

function [hTxt] = fplot_boxplot(subdata, ax, varargin)

fc = [0.6 0.6 0.6;
    mycolors('lightorange');
    mycolors('redorange')];

yl_col = [.65 .65 .65];
yl_width = 2;

boxwidth = 0.1;

labels = {'', '', ''};
if any(strcmp(varargin, 'labels'))
    ind = find(strcmp(varargin, 'labels'));
    labels = varargin{ind+1};
end

xlim = ax.XLim;
p = plot(xlim, [median(subdata(:,1)) median(subdata(:,1))], 'color', yl_col, 'linewidth', yl_width, 'linestyle', ':'); 
p.Color(4) = 0.5;
hold on
boxplot(subdata, 'Notch', 'on',...
    'boxstyle', 'filled', 'colors', fc, 'medianstyle', 'line', ...
    'outliersize', 2, 'symbol', '', 'widths', boxwidth,...
    'labelorientation', 'inline', 'positions', linspace(0.1, 0.3, size(subdata,2)), ...
    'Labels', labels)

% % See all tags of all objects available
% a = get(get(gca,'children'),'children');   % Get the handles of all the objects
% t = get(a{1},'tag')

% loop over the boxes
for i = 1:size(subdata,2)
    
    m = findobj('tag', 'Median');
    m(i).XData = m(i).XData +[0.01 -0.01];
    m(i).Color = 'k';
    m(i).LineWidth = 2;
    
    w = findobj('tag', 'Whisker');
%     w(i).LineStyle = 'none';
    w(i).LineStyle = '-';
    w(i).Color = 'none'; %[.5 .5 .5];

    b = findobj('tag', 'Box');
    b(i).LineWidth = 5;
    
    nh = findobj('tag', 'NotchHi');
    nh(i).Color = nh(i).MarkerEdgeColor;
    nh(i).Marker = 'none';
    nh(i).LineStyle = '-';
    nh(i).LineWidth = 0.5;
    nh(i).XData = nh(i).XData -[0.05 -0.05];
    nh(i).YData = [nh(i).YData nh(i).YData];

    nl = findobj('tag', 'NotchLo');
    nl(i).Color = nl(i).MarkerEdgeColor;
    nl(i).Marker = 'none';
    nl(i).LineStyle = '-';
    nl(i).LineWidth = 0.5;
    nl(i).XData = nl(i).XData -[0.05 -0.05];
    nl(i).YData = [nl(i).YData nl(i).YData];

    uistack(m(i),'top')
    
    % Add semi-transparent dots instead of whiskers
    x = w(i).XData(1);
    idx = fliplr(1:size(subdata,2));
    y = subdata(:,idx(i));
    ybox = b(i).YData;
    to_plot = y>ybox(2) | y<ybox(1);
    s = scatter(repmat(x, 1, sum(to_plot)), y(to_plot), 1.5, nl(i).MarkerEdgeColor);
    s.MarkerFaceAlpha = 0.55;
    s.MarkerEdgeAlpha = 0.55;
    
end


%% Format axis

% set axis to whisker length
q = quantile(subdata, [0.25 0.5 0.75]);
wd = 1.5;
whisker_ranges = [q(1,:) - wd*(q(3,:)-q(1,:)); q(3,:) + wd*(q(3,:)-q(1,:)) ];
max_ranges = [min(min(whisker_ranges)) max(max(whisker_ranges))];
set(gca, 'YLim', max_ranges) % - [0.01*diff(max_ranges) -0.01*diff(max_ranges)])

hTxt=findobj(findobj(gca,'Type','hggroup'),'Type','Text');

set(gca, 'XAxisLocation', 'origin', 'Box', 'off', 'XTickLabels', [], ...
    'XLim', [0.05 0.35])





end