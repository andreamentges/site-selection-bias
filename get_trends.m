%% Gets biodiversity trends
% Returns the detected change as a linear slope estimate

function [detection] = get_trends(samples, varargin)

nsamples = size(samples.S,1);
tend     = size(samples.S,2);

fields = fieldnames(samples);

% Loop over fields (diversity measures)
for f = 1 : length(fields)
    
    field = eval(sprintf('samples.%s', fields{f}));

    eval(sprintf('detection.%s_slope     = NaN(nsamples, 1);', fields{f}));
    eval(sprintf('detection.%s_mean      = NaN(nsamples, 1);', fields{f}));
    eval(sprintf('detection.%s_std       = NaN(nsamples, 1);', fields{f})); 
 
    
        for s = 1 : nsamples
        
            lm = fitlm(1:tend, field(s,:),'linear');
            lmcoeff = table2array(lm.Coefficients);
            eval(sprintf('detection.%s_slope(s) = lmcoeff(2,1);', fields{f}));
            eval(sprintf('detection.%s_mean(s) = mean(field(s,:));', fields{f}));
            eval(sprintf('detection.%s_std(s)  = std(field(s,:));', fields{f}));

            if isinf(lmcoeff(2,3)) && length(unique(field(s,:)))==1 
                eval(sprintf('detection.%s_slope(s)  = 0;', fields{f}));
            end
        end 
 
end
 

end