function Metric = MapMetric(NicheMap, ModelMap, show, R)
% Calculates a metric for the similarity between two maps, NicheMap and ModelMap
% 
% INPUTS:
%   NicheMap: A matrix representing the true niche map
%   ModelMap: A matrix representing the model's predicted niche map
%   show: A flag for whether or not to display the resulting plot
%   R: The geographic reference object of the map (optional)
% 
% OUTPUTS:
%   Metric: A numerical value representing the similarity between the two maps
%
% NOTES:
% - NicheMap and ModelMap must have the same size and dimensions
% - If R is not provided, the resulting plots will not have geographic information
% - The Metric output is a vector that includes both the main metric value and a secondary 1-norm metric value
% - If show is set to 1, a figure will be displayed showing the original and estimated maps, along with their corresponding
%   plots of cumulative sum normalized by the maximum pixel value

    % Check if R is provided, otherwise initialize as empty
    if nargin < 4
        R = [];
    end
    
    % Define NicheMap and ModelMap without NaN values
    VectorNicheMap = NicheMap(:);
    IdxNicheMap = ~isnan(VectorNicheMap);
    VectorNicheMap = VectorNicheMap(IdxNicheMap);
    
    VectorModelMap = ModelMap(:);
    VectorModelMap = VectorModelMap(IdxNicheMap);
    
    % Extract only the true niche map values (pixels where NicheMap and ModelMap are both non-zero)
    combinedMap = VectorNicheMap + VectorModelMap;
    absoluteIdx = combinedMap ~= 0;
    VectorNicheMap = VectorNicheMap(absoluteIdx);
    VectorModelMap = VectorModelMap(absoluteIdx);
    
    % Calculate main metric using 1-norm
    metric2 = 1 - sum(abs(VectorNicheMap-VectorModelMap),'omitnan')/sum(absoluteIdx);
    
    % Sort the model map pixels according to the niche map values
    [VectorNicheMap, IdxModelMap] = sort(VectorNicheMap, 'ascend');
    VectorModelMap = VectorModelMap(IdxModelMap);
    
    % Calculate cumulative sum normalized by the maximum pixel value for each sorted pixel
    VectorModelMap = cumsum(VectorModelMap, 'omitnan');
    VectorNicheMap = cumsum(VectorNicheMap, 'omitnan');
    maxim = max(max(VectorNicheMap), max(VectorModelMap)); 
    VectorNicheMap = VectorNicheMap/maxim;
    VectorModelMap = VectorModelMap/maxim;
    
    % Calculate the difference map
    DiffMap = abs(VectorNicheMap - VectorModelMap);
    LengDiffMap = length(DiffMap);
    LengDiffMap = linspace(0, 1, LengDiffMap);
    
    % Calculate the main metric using trapezoidal integration
    Metric = 1 - trapz(LengDiffMap, DiffMap);
    
    % Display figures if show is set to 1
    if show == 1   
        % Plot niche and model maps
        figure
        LengthNicheMap = length(VectorNicheMap);
        LengthNicheMap = linspace(0, 1, LengthNicheMap); 
        clf
        plot(LengthNicheMap, VectorModelMap, 'LineWidth', 2)
        hold on
        plot(LengthNicheMap, VectorNicheMap, 'LineWidth', 2)
        legend('Estimated','Original','Location','best')
        title(strcat(num2str(round(Metric*100,2)),'%'))

        % Plot niche and model maps on a map
        if ~isempty(R)
            figure
            subplot(1,2,1)
            geoshow(NicheMap, R, 'DisplayType','surface');
            axis off
            title('Niche')
            subplot(1,2,2)
            geoshow(ModelMap, R, 'DisplayType','surface');
            contourcmap('jet', 0 : 0.05 : 1, 'colorbar', 'on', 'location', 'vertical')
            axis off
            title('Model')
        end
    end

    % Append metric2 to Metric vector
    Metric = [Metric, metric2];

end