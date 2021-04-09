%% Make SAD from S and I
% Returns a species abundance distribution of specified type and (optional)
% properties.
%
% Example: get_SAD(Stot, 'lognormal', 'mu', mu, 'var', var)
%
% Note:
% Increasing mu increases the number of individuals, but makes it more
% even.
% Increasing var increases the number of individuals too, but makes it less
% even.


function [Ni] = get_SAD(S, type, varargin)
%% Security checks

assert(rem(S,1)==0, 'S needs to be integer!')
assert(S>0, 'S needs to be >0!')

% Default: no visualization
visualize = 'off';
if any(strcmp(varargin, 'visualize'))
    ind = find(strcmp(varargin, 'visualize'));
    visualize = varargin{ind+1};
end

%% Generate SAD

if strcmp(type, 'uniform') 
    
    % default or specified number of individuals per species
    N = 100;
    if any(strcmp(varargin, 'N'))
        ind = find(strcmp(varargin, 'N'));
        N = varargin{ind+1};
    end
    
    % Get SAD
    Ni = repmat(N, S, 1);

elseif strcmp(type, 'lognormal') 
    
    % mean of lognormal distribution
    mu = 4;
    if any(strcmp(varargin, 'mu'))
        ind = find(strcmp(varargin, 'mu'));
        mu = varargin{ind+1};
    end
     
    % variance of lognormal distribution
    var = 1.5;
    if any(strcmp(varargin, 'var'))
        ind = find(strcmp(varargin, 'var'));
        var = varargin{ind+1};
    end
        
    % get SAD - like in R function rpoilog
    lambda = exp(var * randn([1, S]) + mu + log(1));
    Ni = sort(poissrnd(lambda, [1 S]), 'descend')';
    Ni(Ni==0) = 1;

else
    error('Undefined type for SAD. Valid types: "uniform".')
end

if strcmp(visualize, 'on')
    
    figure('color', 'white', 'position', [-1074,505,314,366])

    % Histogram (individual bins per species)
	subplot(3,1,1)
    edges = 0.5:1:(max(Ni)+0.5);
    h = histogram(Ni, edges, 'edgecolor', mycolors('green'));
    xlabel('N'), ylabel('Count')
    
    % Scatter plot
    subplot(3,1,2)
    Nis = 1:max(Ni);
    plot(Nis(h.Values>0), h.Values(h.Values>0), 'LineStyle', 'none', ...
        'Marker', '.', 'Color', 'r')
    set(gca, 'YLim', [0 inf])
    xlabel('N'), ylabel('Count')
    
    % Rank abundance curve
    subplot(3,1,3)
    yr = sort(Ni/sum(Ni), 'descend');
    plot(yr, 'LineStyle', 'none', 'Marker', '.', 'Color', mycolors('purple'))
    xlabel('Rank'), ylabel(sprintf('Relative\nabundance'))
    set(gca, 'YScale', 'log')
    
    
end

%% Security checks
assert(all(rem(Ni,1)==0), 'All numbers of individuals must be integers!')
assert(all(Ni>0), 'All numbers of individuals must greater zero!')

end
