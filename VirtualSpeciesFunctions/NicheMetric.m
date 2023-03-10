function Metric = NicheMetric(Maps, Dimensions, Method, Plotting)
% Measures the performance of the Virtual Species method
%
% Inputs:
%   maps - A 1-by-k array of map structs, where each struct has a field
%          called 'NormDistance' containing the normalized distance of
%          each point in the map from the nearest niche center.
%   Dimensions - The number of dimensions of the maps.
%   method - The method used to calculate the metric. It can be one of:
%            'ShannonEntropy', 'LorenzCurve', 'KolmogorovSmirnov',
%            'Kurtosis', 'Rank', or 'kkplot'. Default is 'LorenzCurve'.
%   plotting - A logical indicating whether to plot the Lorenz curve or not.
%              Default is true.
%
% Outputs:
%   met - A struct containing the following fields:
%         - Metric: The value of the metric for the given method.
%         - MapMetric: An array of size dim-by-1 containing the metric
%                      values for each dimension.
%         - SortedNormalizedIndex: An array of size dim-by-1 containing the
%                                  metric values for each dimension sorted
%                                  in ascending order.
%         - idx: An array of size dim-by-1 containing the indices of the
%                sorted metric values.
%

if nargin < 3
    Method = 'LorenzCurve';
end

if nargin < 4
    Plotting = true;
end

k = length(Maps);
dim = length(Maps(1).NormDistance);
mapdists = zeros(k, dim);

for i = 1:k
    mapdists(i,:) = Maps(i).NormDistance;
end

switch Method
    case 'ShannonEntropy'
        Entropy = zeros(dim,1);
        for i = 1:dim
            p = mapdists(:,i);
            Entropy(i) = -(1/log2(k)) * sum(p .* log2(p + eps));
        end
    
        unifEntropy = log2(k);
        mapindex = Entropy ./ unifEntropy;
        index = mean(mapindex);

    case 'LorenzCurve'
        Gini = zeros(dim,1);
        for i = 1:dim
            p = mapdists(:,i);
            [f, x] = ecdf(p);
            Gini(i) = 1 - (sum(((x-f).^2))/100) / (sum(((x).^2))/100);
        end
    
        mapindex = Gini;
        index = mean(Gini);
    
    case 'KolmogorovSmirnov'
        KS = zeros(dim,1);
        Hs = zeros(dim,1);
        for i = 1:dim
            p = mapdists(:,i);
            [f, x] = ecdf(p);
            [h, p, ks2stat] = kstest2(f, x);
            KS(i) = 1 - ks2stat;
            Hs(i) = h;
        end
        mapindex = KS;
        index = mean(Hs);
    
    case 'Kurtosis'
        K = zeros(dim,1);
        maxvar = max(var(mapdists,0,2));
        unifKurt = 9/5;
        
        for i = 1:dim
            p = mapdists(:,i);
            varp = var(p);
            kurt = kurtosis(p);
            K(i) = (varp./maxvar)*(abs((kurt-unifKurt)/max([unifKurt,kurt]))+1);
        end
    
        mapindex = K;
        index = mean(K);
    
    case 'Rank'
        Rank = zeros(dim,1);
        Hs = zeros(dim,1);
        Ps = zeros(dim,1);
        
        for i = 1:dim
            p = mapdists(:,i);
            x = rand(length(p),1);
            [ps,h,stats] = ranksum(p,x);
            Rank(i) = stats.ranksum;
            Hs(i) = abs(1-h);
            Ps(i) = 1-ps;
        end
    
        mapindex = Hs;
        index = mean(Hs);
    
    case 'kkplot'
        EK = zeros(dim,1);
        
        for i = 1:dim
            p = mapdists(:,i);
            EK(i) = kurtosis(p)-3;
        end
        
        unifEK = 6/5*ones(dim,1);
        linmod = fitlm(EK, unifEK);
        plot(linmod);
        mapindex = abs(EK./unifEK);
        index = mean(mapindex);
    
end

[SortedIndex,idx] = sort(mapindex);

if Plotting
    Info.idx = idx;
    Info.SortNormDistance = SortedIndex;
    Map = NicheGeneration(Dimensions, Info, 1, true);
end

Metric.Metric = index;
Metric.MapMetric = mapindex;
Metric.SortedNormalizedIndex = SortedIndex;
Metric.idx = idx;

end
