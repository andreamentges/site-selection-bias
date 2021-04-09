%% Generate aggregated spatial distribution of individuals
% 
% Returns the coordinates of individuals used for the initial spatial
% distribution of individuals across the landscape in the "aggregated"
% model variant.
%
% How this works:
% 1. for each species with Ni > 1, draw a number of parents between 1 and
% Ni/10
% 2. for each individual of the species, randomly assign a parent
% 3. randomly distribute parents across the plane (homogeneous poisson
% point process)
% 4. distribute individuals around their parents in a circle with radius r

function [coords_x, coords_y] = get_aggregated_coordinates(SAD, c)

Stot = length(SAD); % Number of species
Itot = sum(SAD); % Number of Individuals

% Radius of aggregates (strongly clumped if small, more random if large)
r = c/20 ;

coords_x = [];
coords_y = [];

%% Loop over species

for s = 1:Stot
    
    Ni = SAD(s);
    
    if Ni>1 % distribute individuals in clumps
        
        % randomly draw the number of parents between 1 and Ni/10
        np = randi(max(1, floor(Ni/10)));
        
        % for each individual, assign a random parent
        ip = randi(np, 1, Ni);
        
        % randomly assign center coordinates to the parents (1:1:100)
        px = randi(c, 1, np);
        py = randi(c, 1, np);
        
        % randomly assign relative polar coordinates to the individuals
        theta = 2*pi*(rand(Ni,1))';
        rho = r*sqrt(rand(Ni,1))';
        
        % convert these relative coordinates from polar to Cartesian
        [ix_rel,iy_rel] = pol2cart(theta,rho);
        
        % for each individual, move its relative cartesian coordinates by 
        % the center of their respective parent and round the result
        ix = round(ix_rel+px(ip));
        iy = round(iy_rel+py(ip));
        
    else % if only one individual exists, randomly place it in the plane
        
        ix = randi(c);
        iy = randi(c);
        
    end
    
    coords_x = [coords_x ix];
    coords_y = [coords_y iy];
    
end

% correct for boundary condition
coords_x = mod(coords_x-1,c)+1;
coords_y = mod(coords_y-1,c)+1;

%% Security checks of result

assert(isequal(size(coords_x), [1 Itot]) && isequal(size(coords_y), [1 Itot]),...
    'Wrong size of coordinate vector')
assert(all([coords_x coords_y]>=1) && all([coords_x coords_y]<=c), ...
    'Boundary condition of coordinates violated! Ouside of cxc range!')

end
