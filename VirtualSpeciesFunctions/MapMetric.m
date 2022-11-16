function Metric = MapMetric(NicheMap, ModelMap, show, R)

    if nargin<4
        R=[];
    end
    
    % Deffining the Niche map without nans
    VectorNicheMap = NicheMap(:);
    IdxNicheMap = ~isnan(VectorNicheMap);
    VectorNicheMap = VectorNicheMap(IdxNicheMap);
    
    % Deffining the Model map without nans
    VectorModelMap = ModelMap(:);
    VectorModelMap = VectorModelMap(IdxNicheMap);
    
    % Extracting the true niche    
    combinedMap = VectorNicheMap + VectorModelMap;
    absoluteIdx = combinedMap ~= 0;
    VectorNicheMap = VectorNicheMap(absoluteIdx);
    VectorModelMap = VectorModelMap(absoluteIdx);
    
    %1-Norm for metric
    metric2 = 1 - sum(abs(VectorNicheMap-VectorModelMap),'omitnan')/sum(absoluteIdx);
    
    % Sort the model map pixels according to the niche map values
    [VectorNicheMap, IdxModelMap] = sort(VectorNicheMap, 'ascend');
    VectorModelMap = VectorModelMap(IdxModelMap);
 
    VectorModelMap = cumsum(VectorModelMap, 'omitnan');
    VectorNicheMap = cumsum(VectorNicheMap, 'omitnan');
    
    % Maximum pixel value between both maps
    maxim = max(max(VectorNicheMap), max(VectorModelMap)); 
    
    % Normalization between 0-1 each sorted pixel
    VectorNicheMap = VectorNicheMap/maxim;
    VectorModelMap = VectorModelMap/maxim;
    
    DiffMap = abs(VectorNicheMap - VectorModelMap);
     
    LengDiffMap = length(DiffMap);
    LengDiffMap = linspace(0, 1, LengDiffMap);
    
    Metric = 1 - trapz(LengDiffMap, DiffMap);
    
    if show == 1   
        figure
        LengthNicheMap = length(VectorNicheMap);
        LengthNicheMap = linspace(0, 1, LengthNicheMap); 
        clf
        plot(LengthNicheMap, VectorModelMap, 'LineWidth', 2)
        hold on
        plot(LengthNicheMap, VectorNicheMap, 'LineWidth', 2)
        legend('Estimated','Original','Location','best')
        title(strcat(num2str(round(Metric*100,2)),'%'))
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
    
    Metric = [Metric, metric2];
    
end
