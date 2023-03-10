function MapInfo = NicheGeneration(Layers, InitialPoint, occupation, plotting)
% NicheGeneration generates a niche map for a virtual species based on
% environmental variables.
%
% INPUTS:
%   Layers: a structure generated by the 'ReadLayers' function.
%   InitialPoint: a structure generated by the 'InitialPoint' function.
%   occupation: an integer with the occupation value ([0,1)).
%   show: boolean variable (true, false) to show the resulting niche map.
%
% OUTPUTS:
%   MapInfo: a structure with the following fields:
%       - Map: a matrix representing the niche map.
%       - SortNormDistance: a vector of sorted and normalized distances.
%       - NormDistance: a vector of normalized distances.

    % Start timer
    tic 

    % Read climate variables
    indicator = Layers.Indicator;  % Indicator matrix
    dimension = Layers.Dimensions(1);  % Dimension of the map
    r = Layers.R;  % Geographic reference object
    [Dim1, Dim2] = size(Layers.Map);  % Size of the map
    Map = nan(Dim1, Dim2);  % Initialize the niche map

    % Extract initial point information
    idx = InitialPoint.idx;  % Index of the initial point
    SortNormDistance = InitialPoint.SortNormDistance;  % Sorted and normalized distances

    % Calculate the limit for the niche map
    limit = round(dimension * occupation); 

    % Normalize the distances and set distances beyond the limit to 0
    NormDistance = zeros(1,length(SortNormDistance));
    SortNormDistance(limit : end) = 0;
    IndexSortNorm = 1 : limit - 1;
    SortNormDistance(IndexSortNorm) = normalize(SortNormDistance(IndexSortNorm), 2, 'range');
    
    % Assign the normalized distances to the niche map
    NormDistance(idx) = SortNormDistance;
    Map(~indicator) = NormDistance;

    % Plot niche map if show is true
    if plotting
        clf
        geoshow(Map, r, 'DisplayType', 'surface');
        contourcmap('jet', 0 : 0.05 : 1, 'colorbar', 'on', 'location', 'vertical')
    end

    % Save output in a structure
    MapInfo.Map = Map;
    MapInfo.SortNormDistance = SortNormDistance;
    MapInfo.NormDistance = NormDistance;
    
    % Stop timer
    toc      
end
