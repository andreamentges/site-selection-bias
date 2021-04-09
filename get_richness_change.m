%% Change species identities to change richness without changing N
% Imposes the biodiversity change regime "richness change" for the
% specified list of species, specified sampling duration and specified rate
% of biodiversity change.

function [species] = get_richness_change(species_in, tend, change_rate)

% total number of individuals (is constant)
I = size(species_in,1);

% Species-identity of individuals over time matrix
species = NaN(I, tend);
species(:,1) = species_in;

% check whether change rate can be realized
% % % if change_rate>0
    Ni = hist(species_in, unique(species_in));
    available_slots = Ni-1;
    rate_ok = sum(available_slots)>=abs(change_rate)*tend;
    assert(rate_ok, 'Richness change cannot be realized. Increase number of individuals, decrease change_rate, or decrease tend.')

    % "change_rate" individuals of the most common species are turned into new species
    for t = 2:tend
        Ni = hist(species(:,t-1), unique(species(:,t-1)));
        available_slots = Ni-1;
        [~, idx] = sort(available_slots, 'ascend');
        newspecs = max(species(:,t-1)) + [1:abs(change_rate)];
        losers   = fliplr(idx);
        species(:,t) = species(:, t-1);
        placed = 0;
        left = abs(change_rate);
        i = 1;
        while left>0
            nplacable = min(left, available_slots(losers(i)));
            species(find(species(:,t)==losers(i), nplacable , 'first'),t) = newspecs(1:nplacable);
            placed = placed + nplacable;
            left = abs(change_rate)-placed;
            newspecs(1:nplacable) = [];
            i = i+1;
        end    
    end
    
% for negative change rates: just flip increase time series
if change_rate < 0
    species = fliplr(species);
end

assert(all(sum(species>0)==sum(species(:,1)>0)), 'S changes')

