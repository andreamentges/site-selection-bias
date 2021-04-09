%% Construct changing species identities
% Imposes the biodiversity change regime "evenness change" for the
% specified list of species, specified sampling duration and specified rate
% of biodiversity change.

function [species] = get_evenness_change(species_in, tend, change_rate)

I = size(species_in,1);

species = NaN(I, tend);
species(:,1) = species_in;

% mean number of individuals per species
Imean = length(species_in)/max(species_in);

if change_rate > 0 % make more even
    % "change_rate" individuals of the most common species are turned into the rarest species
    for t = 2:tend
        species(:,t) = species(:, t-1);
        Ni = hist(species(:,t), unique(species(:,t)));
        Idiff = Ni-Imean; % difference to the mean (perfectly even)
        [~, diffidx] = sort(Idiff); % first entry: rarest species
        [~, winners] = find(Idiff==min(Idiff));  
        winner = randsample(winners, change_rate, 'true');
        placed = 0;
        left = change_rate;
        i = 0;
        while left>0
            loser = diffidx(end-i);
            n2place = min(Ni(loser)-1, left);
            species(find(species(:,t)==loser, n2place , 'first'), t) = winner(1:n2place);
            placed = placed + n2place;
            left = change_rate-placed;
            i = i+1; 
        end
    end
elseif change_rate < 0 % make less even
    % "change_rate" individuals of the next-rarest species are turned into the most common species
    warning('Evenness change is extremely slow at the moment! UPDATE!')

    for t = 2:tend
        species(:,t) = species(:, t-1);
        for i = 1:abs(change_rate)
            Ni = hist(species(:,t), unique(species(:,t)));
            losers = find(Ni == min(Ni(Ni>1)));
            [~, winner] = max(Ni);
            species(find(species(:,t)==losers(randperm(length(losers),1)), 1, 'first'),t) = winner;
        end        
    end
        
else
    error('change rate cannot be 0 for evenness change')
end

assert(all(sum(species>0)==sum(species(:,1)>0)), 'S changes')
