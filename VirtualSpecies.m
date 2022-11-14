function MapInfo = VirtualSpecies(ReadInfo, InfoInitialPoint, Occupation, show)
        
    %%%-- Reading clim variables--%%

    Indicator = ReadInfo.Indicator;
    Dimension = ReadInfo.Dimensions(1);
    R = ReadInfo.R;
    [Dim1, Dim2] = size(ReadInfo.Map);
    Map = nan(Dim1, Dim2);
        
    NormDistance = InfoInitialPoint.NormDistance;
    idx = InfoInitialPoint.idx;
    SortNormDistance = InfoInitialPoint.SortNormDistance;
    
    limit = round(Dimension * Occupation); 
    
    SortNormDistance(limit : end) = 0;
    IndexSortNorm = 1 : limit - 1;
    SortNormDistance(IndexSortNorm) = normalize(SortNormDistance(IndexSortNorm), 2, 'range');
    
    NormDistance(idx) = SortNormDistance;
    Map(~Indicator) = NormDistance;
    
    if show == 1
        figure(1)
        clf
        geoshow(Map, R, 'DisplayType', 'surface');
        contourcmap('jet', 0 : 0.05 : 1, 'colorbar', 'on', 'location', 'vertical')
    end
    
    MapInfo = InfoInitialPoint;
    MapInfo.Map = Map;
    MapInfo.SortNormDistance = SortNormDistance;
    MapInfo.NormDistance = NormDistance;
      
end
