%% Make observation in virtual landscape
%
% Returns the diversity timeseries observed within a selected sampling
% grain.

function [div, grain2mean] = get_diversity_timeseries(Stot, c, g, tend, nsamples, change, change_rate_rel, varargin)
%% Preparations & checks

assert(isequal(size(Stot), [1 1]), 'Specify Stot as first input argument.')

% Default: no visualization of sampled grain
visualize = 'off';
if any(strcmp(varargin, 'visualize')) || any(strcmp(varargin, 'plot_grain'))
    ind = find(strcmp(varargin, 'visualize'));
    visualize = varargin{ind+1};
end

% Default: no plotting of movement
movement = 'off';
if any(strcmp(varargin, 'movement')) || any(strcmp(varargin, 'plot_movement'))
    ind = find(strcmp(varargin, 'movement'));
    movement = 'on';
end

% Default: no plotting of movement
plot_diversity = 'off';
if any(strcmp(varargin, 'plot_diversity'))
    plot_diversity = 'on';
end

% Default: random sampling. Alternatively: site-selection bias.
sampling = 'random';
if any(strcmp(varargin, 'sampling'))
    ind = find(strcmp(varargin, 'sampling'));
    sampling = varargin{ind+1};
end

% Default: no timelag, all time points used
timelag = 0;
if any(strcmp(varargin, 'timelag'))
    ind = find(strcmp(varargin, 'timelag'));
    timelag = varargin{ind+1};
end
if timelag == 0
    tendL = tend;
else
    tendL = tend*timelag;
end

% Default: no truncating. Alternatively: truncate pt first time points of timeseries
truncate   = 'off';
pt = 0;
if any(strcmp(varargin, 'truncate'))
    ind = find(strcmp(varargin, 'truncate'));
    pt = varargin{ind+1};
    if pt == 0
        truncate = 'off';
    else
        truncate = 'on';
        pt = varargin{ind+1}; 
    end
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

% initial spatial distribution of individuals
start_coordinates = 'random';
if any(strcmp(varargin, 'aggregated'))
    start_coordinates = 'aggregated';
elseif any(strcmp(varargin, 'heterogeneous'))
    start_coordinates = 'heterogeneous';    
end
 

%% Pre-assign empty vectors

% These include only observed time points (minus timelag)
div.S     = NaN(nsamples, tend);
div.N     = NaN(nsamples, tend);
div.S_PIE = NaN(nsamples, tend);

for s = 1 : nsamples
    %% Get a random SAD

    % Get species-abundance-distribution (Ni)
    % Check whether the variables mu and var exist
    if exist('mu', 'var') && exist('var', 'var')
        SAD = get_SAD(Stot, 'lognormal', 'mu', mu, 'var', var, 'visualize', 'off');
    elseif ~exist('var', 'var') && ~exist('var', 'var')
        SAD = get_SAD(Stot, 'lognormal', 'visualize', 'off');
    else
        error('Did you specify only one of mu and var for SAD?')
    end

    % for each individual, give species as integer number
    species_in = repelem(1:Stot, SAD)';
    
    %% Simulate movement

    % get absolute change rate (number of individuals) from relative change
    % rate
    change_rate = round(change_rate_rel*sum(SAD));
    
    % Simulation
    if strcmp(start_coordinates, 'heterogeneous') 
        % In this case, the SAD changes because it is constructed from 4
        % individual SADs
        [posXL, posYL, SAD] = simulate_timeseries(SAD, c, g, tendL, change, change_rate, start_coordinates);
        species_in = repelem(1:Stot, SAD)';
    else
        [posXL, posYL] = simulate_timeseries(SAD, c, g, tendL, change, change_rate, start_coordinates);
    end
    
    % In the case of "evenness change", the basic movement is default, but
    % the change is calculated separately as a changing SAD
    if strcmp(change, 'evenness')
        species = get_evenness_change(species_in, tend, change_rate);
    elseif strcmp(change, 'S')
        species = get_richness_change(species_in, tend, change_rate);
    else
        species = repmat(species_in, 1, tend);
    end
    
    %% Get sampling grain

    % find lower-left point of all possible gxg rectangles within cxc
    % (lower left points range between 0 and c-g)
    num_grains = (floor(c/g)^2);
    gsteps = (0:((c/g)-1))*g;
    [gx_all_mesh, gy_all_mesh] = meshgrid(gsteps, gsteps);
    gx_all = gx_all_mesh(:);
    gy_all = gy_all_mesh(:);
    
    % get the initial abundance of all grains
    comm_t0_N = zeros(num_grains,1);
    comm_t0_x = posXL(:,1);
    comm_t0_y = posYL(:,1);
    for i = 1:length(comm_t0_x)
        ind = (gx_all == floor(comm_t0_x(i)/g)*g) & ...
        (gy_all == floor(comm_t0_y(i)/g)*g);
        comm_t0_N(ind) = comm_t0_N(ind) + 1;
    end
    
    % get the initial richness of all grains
    visual = get_visual_specifications(species);
    rich_t0 = zeros(num_grains,1); % richness per grain
    for i = 1:Stot
        occs = unique([floor(posXL(species(:,1)==i,1)/g)*g ...
            floor(posYL(species(:,1)==i,1)/g)*g], 'rows');
        ind = ismember([gx_all gy_all], occs, 'rows');
        rich_t0(ind) = rich_t0(ind) + 1;
    end

    % Choose sampling grain according to sampling strategy
    if strcmp(sampling, 'random') % Random grain
        rand_idx = randi([1 num_grains], 1, 1);
        gx = gx_all(rand_idx);
        gy = gy_all(rand_idx);
        
    elseif strcmp(sampling, 'pop-biased') || strcmp(sampling, 'biased') % Population bias
        % find the rectangle among these with the highest initial abundance 
        % of species 1
        spec1_t0_x = posXL(species(:,1)==1,1);
        spec1_t0_y = posYL(species(:,1)==1,1);
        spec1_t0_N = zeros(num_grains,1);
        for i = 1:length(spec1_t0_x)
            ind = (gx_all == floor(spec1_t0_x(i)/g)*g) & ...
                (gy_all == floor(spec1_t0_y(i)/g)*g);
            spec1_t0_N(ind) = spec1_t0_N(ind) + 1;
        end
        [~,idx] = sort(spec1_t0_N, 'descend');
        gx = gx_all(idx(1));
        gy = gy_all(idx(1));
       
    elseif strcmp(sampling, 'comm-biased') % Community bias
        % find the rectangle with the highest initial abundance 
        % of species (all species, total N in the grain)
        [~,idx] = sort(comm_t0_N, 'descend');
        gx = gx_all(idx(1));
        gy = gy_all(idx(1));   
    elseif strcmp(sampling, 'rich-biased') % Richness bias
        % find the rectangle among these with the highest initial S 
        [~,idx] = sort(rich_t0, 'descend');
        gx = gx_all(idx(1));
        gy = gy_all(idx(1)); 
    else
        error('Unspecified sampling strategy. Choose "random" or "biased".')
    end
    
    %% Save in which quadrant the chosen grain lies
    
    grain2mean.chosen_quadrant(s) = 1*(gx<c/2 && gy>=c/2) ...
        + 2*(gx>=c/2 && gy>=c/2) ...
        + 3*(gx<c/2 && gy<c/2) ...
        + 4*(gx>=c/2 && gy<c/2);

    %% Visualize the sampled grain
    if strcmp(visualize, 'on')
        visual = get_visual_specifications(species);
        figure('color', 'white', 'position', [-548,464,342,276])
            hold on
            set(gca, 'XLim', [0 c], 'YLim', [0 c])
        xlabel('Horizontal coordinate')
        ylabel('Vertical coordinate')
        if strcmp(sampling, 'pop-biased') || strcmp(sampling, 'biased')
            splot = scatter(posXL(:,1), posYL(:,1), visual.szs(:,1), visual.colorder(species(:,1),:), 'filled');
            splot.MarkerFaceAlpha = 0.5; 
            scatter(posXL(species(:,1)==1,1), posYL(species(:,1)==1,1), 'r')
            title('Sampling grain (red: species 1 at t=0)')
        elseif strcmp(sampling, 'random') || strcmp(sampling, 'comm-biased') ...
                 || strcmp(sampling, 'rich-biased')
            scatter(posXL(:,1), posYL(:,1), visual.szs(:,1), visual.colorder(species(:,1),:), 'filled');
        end
        plot(gx+[0 g g 0 0], gy+[0 0 g g 0], 'r', 'linewidth', 1.5)
        set(gca, 'XTick', 0:g:c, 'YTick', 0:g:c)
        grid on
        % Cannot change the xticks when grid is on. Plot your own grid!
        
    end

    %% Make virtual observation
    
    % If "timelag" is specified, reduce time series respectively (so length
    % of time series corresponds to tend in the end)
    if timelag > 0 
        keep = 1:timelag:tendL;
    else
        keep = 1:tendL;
    end
    posX = posXL(:,keep);
    posY = posYL(:,keep);

    % Pre-assign vectors for changing N, Ni,.. in every time step
    Stot_max = max(max(species));
    Ni    = NaN(Stot_max,tend);
    N     = NaN(1,tend);
    pi    = NaN(Stot_max,tend);
    S     = NaN(1,tend);
    PIE   = NaN(1,tend);
    S_PIE = NaN(1,tend);

    % Boolean: is the individual within or without grain?
    isg = posX>=gx & posX<gx+g & posY>=gy & posY<gy+g;

    % Diversity metrics over time
    for t = 1 : tend
        Ni(1:max(species(:,t)),t) = accumarray(species(:,t), isg(:,t));
        Ni(isnan(Ni(:,t)),t) = 0;
        N(t)     = sum(Ni(:,t));
        pi(:,t)  = Ni(:,t)/N(t);
        S(t)     = sum(Ni(:,t)>0);
        PIE(t)   = 1-sum(pi(:,t).^2);
        S_PIE(t) = 1/(1-PIE(t));  
    end   
    

    %%%%%%%%%
    % Get true diversity in landscape over time!
    %%%%%%%%%
    % Boolean: is the individual present at all (or NaN because it will be
    % added later due to change)
    isp = ~isnan(posX);
    % Pre-assign vectors
    tNi    = NaN(Stot_max,tend);
    tN     = NaN(1,tend);
    tpi    = NaN(Stot_max,tend);
    tS     = NaN(1,tend);
    tPIE   = NaN(1,tend);
    tS_PIE = NaN(1,tend);
    
    % Diversity metrics over time
    for t = 1 : tend
        tNi(1:max(species(:,t)),t)  = accumarray(species(:,t), isp(:,t));
        tNi(isnan(tNi(:,t)),t) = 0;
        tN(t)     = sum(tNi(:,t));
        tpi(:,t)  = tNi(:,t)/tN(t);
        tS(t)     = sum(tNi(:,t)>0);
        tPIE(t)   = 1-sum(tpi(:,t).^2);
        tS_PIE(t) = 1/(1-tPIE(t));  
    end 
    
    % plot movement for comparison
    if strcmp(movement, 'on')
        % Plot the movement with:
        plot_movement(SAD, c, g, tend, posX, posY, species, 'wait');
    end

    div.N(s,:)     = N;
    div.S(s,:)     = S;
    div.S_PIE(s,:) = S_PIE;
    div.tN(s,:)     = tN;
    div.tS(s,:)     = tS;
    div.tS_PIE(s,:) = tS_PIE;

if strcmp(plot_diversity, 'on')
    plot_diversity_dynamics(div)
end

    %% Get spatial heterogeneity:
    % initial (at time = 1) difference of S (or N) in the selected site to 
    % mean S (or N) in all other sites

    % difference to mean abundance
    grain2mean.N(s,1) = N(1)-mean(comm_t0_N);
    grain2mean.N_st(s,1) = (N(1)-mean(comm_t0_N))/sum(SAD);
    % assert(sum(SAD)/num_grains==mean(comm_t0_N), 'not the same, should be the same!')
    
    % difference to mean richness
    grain2mean.S(s,1) = S(1)-mean(rich_t0);
    grain2mean.S_st(s,1) = (S(1)-mean(rich_t0))/Stot;


end

% Optionally left-truncate the time series (cut first pt percent of
% points)
if strcmp(truncate, 'on')  

    div.N(:,1:pt)      = [];
    div.S(:,1:pt)      = [];
    div.S_PIE(:,1:pt)  = [];       
    div.tN(:,1:pt)     = [];
    div.tS(:,1:pt)     = [];
    div.tS_PIE(:,1:pt) = [];

end
    

end