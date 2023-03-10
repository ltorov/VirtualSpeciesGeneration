function Metric = MapMetric(niche_map, model_map, plotting, geo_ref)
% Calculates a metric for the similarity between two maps, NicheMap and ModelMap
% 
% INPUTS:
%   niche_map: A matrix representing the true niche map
%   model_map: A matrix representing the model's predicted niche map
%   plotting: A flag for whether or not to display the resulting plot
%   geo_ref: The geographic reference object of the map (optional)
% 
% OUTPUTS:
%   Metric: A numerical value representing the similarity between the two maps
%
% NOTES:
% - niche_map and model_map must have the same size and dimensions
% - If geo_ref is not provided, the resulting plots will not have geographic information
% - The Metric output is a vector that includes both the main metric value and a secondary 1-norm metric value
% - If plotting is set to 1, a figure will be displayed showing the original and estimated maps, along with their corresponding
%   plots of cumulative sum normalized by the maximum pixel value

    % Check if R is provided, otherwise initialize as empty
    if nargin < 4
        geo_ref = [];
    end
    
    % Define niche_map and model_map without NaN values
    niche_map_vector = niche_map(~isnan(niche_map));
    model_map_vector = model_map(~isnan(niche_map));
        
    % Extract only the true niche map values (pixels where niche_map and model_map are both non-zero)
    combined_map = niche_map_vector + model_map_vector;
    absolute_idx = combined_map ~= 0;
    niche_map_vector = niche_map_vector(absolute_idx);
    model_map_vector = model_map_vector(absolute_idx);

    % Calculate main metric using 1-norm
    metric2 = 1 - sum(abs(niche_map_vector - model_map_vector),'omitnan')/sum(absolute_idx);
    
    % Sort the model map pixels according to the niche map values
    [niche_map_vector, model_map_idx] = sort(niche_map_vector, 'ascend');
    model_map_vector = model_map_vector(model_map_idx);
    
    % Calculate cumulative sum normalized by the maximum pixel value for each sorted pixel
    model_map_vector = cumsum(model_map_vector, 'omitnan');
    niche_map_vector = cumsum(niche_map_vector, 'omitnan');
    maxim = max(max(niche_map_vector), max(model_map_vector)); 
    niche_map_vector = niche_map_vector/maxim;
    model_map_vector = model_map_vector/maxim;
    
    % Calculate the difference map
    diff_map = abs(niche_map_vector - model_map_vector);
    diff_map_len = linspace(0, 1, length(diff_map));
    
    % Calculate the main metric using trapezoidal integration
    metric = 1 - trapz(diff_map_len, diff_map);
    
    % Display figures if show is set to 1
    if plotting == 1   
        % Plot niche and model maps
        figure
        niche_map_len = linspace(0, 1, length(niche_map_vector)); 
        clf
        plot(niche_map_len, model_map_vector, 'LineWidth', 2)
        hold on
        plot(niche_map_len, niche_map_vector, 'LineWidth', 2)
        legend('Estimated','Original','Location','best')
        title(strcat(num2str(round(metric*100,2)),'%'))

        % Plot niche and model maps on a map
        if ~isempty(geo_ref)
            figure
            subplot(1,2,1)
            geoshow(niche_map, geo_ref, 'DisplayType','surface');
            axis off
            title('Niche')
            subplot(1,2,2)
            geoshow(model_map, geo_ref, 'DisplayType','surface');
            contourcmap('jet', 0 : 0.05 : 1, 'colorbar', 'on', 'location', 'vertical')
            axis off
            title('Model')
        end
    end

    % Append metric2 to Metric vector
    Metric = [metric, metric2];

end