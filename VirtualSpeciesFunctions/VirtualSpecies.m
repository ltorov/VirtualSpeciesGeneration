function MapInfo = VirtualSpecies(ReadInfo, InfoInitialPoint, Occupation, show)
    % MapInfo = VirtualSpecies(ReadInfo, InfoInitialPoint, Occupation, show)
    %
    % DESCRIPTION:
    % VirtualSpecies generates a virtual species distribution map.
    % The function reads the climate variables from ReadInfo, and uses
    % the initial point from InfoInitialPoint to generate the map.
    % The function then returns a struct containing the map, as well as
    % updated version of InfoInitialPoint.
    %
    % INPUTS:
    %   ReadInfo: A struct containing information about the climate
    %             variables used to generate the map.
    %   InfoInitialPoint: A struct containing the initial point used to
    %                     generate the map.
    %   Occupation: A scalar value between 0 and 1 specifying the amount
    %               of the climate variable to be used.
    %   show: A boolean value (true/false) that indicates whether or not to
    %         display the generated map.
    %
    % OUTPUTS:
    %   MapInfo: A struct containing the updated version of
    %            InfoInitialPoint, as well as the generated map.
    
    % Read in clim variables
    Indicator = ReadInfo.Indicator;
    Dimension = ReadInfo.Dimensions(1);
    R = ReadInfo.R;
    [Dim1, Dim2] = size(ReadInfo.Map);
    Map = nan(Dim1, Dim2);
    
    % Set up the virtual species distribution map
    NormDistance = InfoInitialPoint.NormDistance;
    idx = InfoInitialPoint.idx;
    SortNormDistance = InfoInitialPoint.SortNormDistance;
    
    % Normalize the top "Occupation" percent of values to 1
    limit = round(Dimension * Occupation);
    SortNormDistance(limit : end) = 0;
    IndexSortNorm = 1 : limit - 1;
    SortNormDistance(IndexSortNorm) = normalize(SortNormDistance(IndexSortNorm), 2, 'range');
    
    % Update NormDistance and Map
    NormDistance(idx) = SortNormDistance;
    Map(~Indicator) = NormDistance;
    
    % Display the map if show is true
    if show
        figure(1)
        clf
        geoshow(Map, R, 'DisplayType', 'surface');
        contourcmap('jet', 0 : 0.05 : 1, 'colorbar', 'on', 'location', 'vertical')
    end
    
    % Output struct
    MapInfo = InfoInitialPoint;
    MapInfo.Map = Map;
    MapInfo.SortNormDistance = SortNormDistance;
    MapInfo.NormDistance = NormDistance;
end
