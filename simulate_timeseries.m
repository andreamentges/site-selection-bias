%% Simulate individuals moving trough landscape over time and in space
% in each time step, change_rate individuals are added, until the input
% SAD is reached at t=tend.
%
% S         total number of species in the landscape
% I         total number of individuals in the landscape
% c         number of cells in grid (landscape size)
% g         number of cells in sampling grain (grain size)
% change_type = 'none', 'N', 'S', or 'evenness'
% change_rate = individuals removed/added per time point
%
% Optional input arguments:
% - 'aggregated', for aggregated spatial distribution of individuals at
% start

function [x, y, Ni] = simulate_timeseries(Ni, c, g, tend, change_type, change_rate, varargin)
 

if strcmp(change_type, 'none') && change_rate~=0
    warning('As change type is "none", no change will be implemented in spite of change rate specified.')
    change_rate = 0;
end

%% Security checks

assert(rem(change_rate,1)==0, 'Change rate needs to be integer!') 
assert(rem(c,1)==0 & rem(g,1)==0, 'c and g need to be integers!')
assert(c>=g, 'c needs to be larger than (or at least equal to) g!')
 
%% Define parameters (as default, or change if specified)
 
% probability of moving to random neighbouring field
p_move = 8/9; %0.8;
if any(strcmp(varargin, 'p'))
    ind = find(strcmp(varargin, 'p'));
    p_move = varargin{ind+1};
end

% initial spatial distribution of individuals
start_coordinates = 'random';
if any(strcmp(varargin, 'aggregated'))
    start_coordinates = 'aggregated';
elseif any(strcmp(varargin, 'heterogeneous'))
    start_coordinates = 'heterogeneous';
end

%% If start_coordinates are heterogeneous, generate 4 separate SADs

if strcmp(start_coordinates, 'heterogeneous')
    SAD1 = get_SAD(100, 'lognormal',  'mu', 3.5, 'var', 1.5); 
    SAD2 = get_SAD(100, 'lognormal',  'mu', 1.2, 'var', 1.5); 
    SAD3 = get_SAD(50, 'lognormal',  'mu', 3.5, 'var', 1.5); 
    SAD4 = get_SAD(50, 'lognormal',  'mu', 1.2, 'var', 1.5); 
    Ni = SAD1 + SAD2 + [SAD3; zeros(50,1)] + [SAD4; zeros(50,1)];
    [sorted_Ni] = sort(Ni, 'descend');
    assert(isequal(Ni, sorted_Ni), 'error in Ni sorting')
end

%% Richness change types

if strcmp(change_type, 'none') && change_rate==0 || strcmp(change_type, 'N') && change_rate>=0
    %%%%%%%%%%%%%%%%%
    % NO CHANGE, or INCREASE N
    %%%%%%%%%%%%%%%%%
    Itot_end = sum(Ni);
    ninc = change_rate*(tend-1);
    Itot_start = Itot_end - ninc;
    assert(Itot_start>0, 'Rate of change times time steps exceeds total number of individuals')

elseif strcmp(change_type, 'N') && change_rate<0
    %%%%%%%%%%%%%%
    % DECREASE N
    %%%%%%%%%%%%%%
    Itot_start = sum(Ni);
    ninc = change_rate*(tend-1);
    Itot_end = Itot_start - abs(ninc);
    assert(Itot_end>0, 'Rate of change times time steps exceeds total number of individuals')
    
elseif strcmp(change_type, 'evenness') || strcmp(change_type, 'S')
    %%%%%%%%%%%%%%
    % No adding of individuals: evenness or richness change
    %%%%%%%%%%%%%%
    % Note: the swapping of individual-to-species-assignments takes place
    % later in the get_diversity_timeseries using the function
    % get_evenness_change.m and get_richness_change.m. 
    % Here, no change of positions is implemented,
    % as indivduals remain where they are but only their identities are
    % changed.
    change_rate = 0;
    Itot_end = sum(Ni);
    ninc = change_rate*(tend-1);
    Itot_start = Itot_end - ninc;
    assert(Itot_start>0, 'Rate of change times time steps exceeds total number of individuals')
    
else
    error('Type of change "%s" not defined.', change_type)
end
  

 
%% randomly assign individuals across grid cells
 
% Pre-assign matrices for simulation results
x = NaN(max(Itot_end, Itot_start), tend); % x coordinate
y = NaN(max(Itot_end, Itot_start), tend); % y coordinate 
% construct identity index vector for individuals which are present at
% start (randomly distribute among rows)
all_inds   = randperm(max(Itot_end, Itot_start));
start_inds = all_inds(1:Itot_start);

 
switch start_coordinates
    case 'random'
        % random starting coordinate for each individual of each species
        % !! Note that positions are 1:1:100, while at the end of this function, 
        % valid positions are transformed to [1:1:c]-0.5
        start_x = ceil(rand(Itot_start,1)*(c)); 
        start_y = ceil(rand(Itot_start,1)*(c));
    case 'aggregated'
        % aggregated starting coordinate, first get SAD of starting individuals
        species_all = repelem(1:length(Ni), Ni)';
        species_start = species_all(sort(start_inds, 'ascend'));
        [SAD_start] = hist(species_start,unique(species_start));
        [start_x, start_y] = get_aggregated_coordinates(SAD_start, c);
    case 'heterogeneous'
        % random starting coordinate for each individual of the four SADs
        % in quadrants of the grid
        SAD1_x = ceil(rand(sum(SAD1),1)*(c/2));
        SAD1_y = (c/2) + ceil(rand(sum(SAD1),1)*(c/2));
        SAD2_x = (c/2) + ceil(rand(sum(SAD2),1)*(c/2));
        SAD2_y = (c/2) + ceil(rand(sum(SAD2),1)*(c/2));
        SAD3_x = ceil(rand(sum(SAD3),1)*(c/2));
        SAD3_y = ceil(rand(sum(SAD3),1)*(c/2));
        SAD4_x = (c/2) + ceil(rand(sum(SAD4),1)*(c/2));
        SAD4_y = ceil(rand(sum(SAD4),1)*(c/2));
        start_x = [SAD1_x; SAD2_x; SAD3_x; SAD4_x] ; 
        start_y = [SAD1_y; SAD2_y; SAD3_y; SAD4_y] ;
        
        % remove some starting coordinates, only leave starting individuals
        start_x(randperm(length(start_x), length(start_x)-Itot_start)) = [];
        start_y(randperm(length(start_y), length(start_y)-Itot_start)) = [];

end
 
% fill positions of starting individuals into matrix
x(start_inds,1) = start_x;
y(start_inds,1) = start_y;
 
if change_rate == 0 || (strcmp(change_type, 'N') && change_rate > 0) || strcmp(change_type, 'evenness') || strcmp(change_type, 'S') 
    
    %%%%%%%%%%%%%%%
    % NO CHANGE or INCREASE N or Evenness change or Richness change
    %%%%%%%%%%%%%%%
    
    % construct index vector of individuals which will be placed due to
    % increase (leave first column empty, because no change in t=1)
    inc_inds = all_inds(Itot_start+1:end);
    inc_idx = [NaN(change_rate, 1) reshape(inc_inds, [change_rate, tend-1])];

    % random starting coordinate for the newly placed individuals
    inc_x = [NaN(change_rate, 1) ceil(rand(change_rate,tend-1)*(c))];
    inc_y = [NaN(change_rate, 1) ceil(rand(change_rate,tend-1)*(c))];

    % Run the simulation (let individuals move)
    for t = 2:tend

        % draw random numbers for deciding whether individuals move at all, 
        % right or left, up or down
        rdecide = rand(Itot_end,3);

        % randomly decide whether an individual moves or not
        moves = rdecide(:,1)<=p_move;

        % randomly decide to which direction
        moves_right = rdecide(:,2)<0.5;
        moves_up    = rdecide(:,3)<0.5;

        % get new x and y coordinates
        new_x = x(:,t-1) + moves_right.*moves - ~moves_right .*moves;
        new_y = y(:,t-1) + moves_up.*moves - ~moves_up .*moves;

        % correct for periodic boundary conditions
        x(:,t) = mod(new_x-1,c)+1;
        y(:,t) = mod(new_y-1,c)+1;

        % insert new individuals to random places
        x(inc_idx(:,t),t) = inc_x(:,t);
        y(inc_idx(:,t),t) = inc_y(:,t);

    end
    
    
elseif strcmp(change_type, 'N') && change_rate < 0
    
    %%%%%%%%%%%%%%%
    % DECREASE N
    %%%%%%%%%%%%%%%
    
    % construct index vector of individuals which will be removed (leave 
    % first column empty, because no change in t=1)
    dec_inds = all_inds(Itot_end+1:end);
    dec_idx = [NaN(abs(change_rate), 1) reshape(dec_inds, [abs(change_rate), tend-1])];

    % Run the simulation (let individuals move)
    for t = 2:tend

        % draw random numbers for deciding whether individuals move at all, 
        % right or left, up or down
        rdecide = rand(Itot_start,3);

        % randomly decide whether an individual moves or not
        moves = rdecide(:,1)<=p_move;

        % randomly decide to which direction
        moves_right = rdecide(:,2)<0.5;
        moves_up    = rdecide(:,3)<0.5;

        % get new x and y coordinates
        new_x = x(:,t-1) + moves_right.*moves - ~moves_right .*moves;
        new_y = y(:,t-1) + moves_up.*moves - ~moves_up .*moves;

        % correct for periodic boundary conditions
        x(:,t) = mod(new_x-1,c)+1;
        y(:,t) = mod(new_y-1,c)+1;

        % insert new individuals to random places
        x(dec_idx(:,t),t) = NaN;
        y(dec_idx(:,t),t) = NaN;

    end
    

else
    
    error('not implemented yet.')
    
end
 
% Valid positions are [1:1:c]-0.5
x = x-0.5;
y = y-0.5;
 
end
 
 

