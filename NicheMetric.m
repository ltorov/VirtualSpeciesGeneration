function met = metrica(maps,Dimensions, method, plotting)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
if nargin<3
    method = 'LorenzCurve';
end

if nargin <4
    plotting = true;
end
% method = 'ShannonEntropy','LorenzCurve','KolmogorovSmirnov','Rank','kkplot'
k = length(maps);


dim = length(maps(1).NormDistance);

mapdists = zeros(k,dim);

for i = 1:k
    mapdists(i,:) = maps(i).NormDistance;
end


%map1dist = map1.NormDistance; map2dist = map2.NormDistance;
%map3dist = map3.NormDistance;


%mapdists = [map1dist; map2dist; map3dist];



if strcmp(method,'ShannonEntropy')

    Entropy = zeros(dim,1);
    for i = 1:dim
        p = mapdists(:,i);
        Entropy(i) = -(1/log2(k))*sum(p.*log2(p+eps));
    end

    
    unifEntropy = log2(k);
    mapindex = Entropy./unifEntropy;
    index = mean(mapindex);

end
if strcmp(method,'LorenzCurve')

    Gini = zeros(dim,1);
    for i = 1:dim
        p = mapdists(:,i);
        [f,x] = ecdf(p);
        Gini(i) = 1-(sum(((x-f).^2))/100)/(sum(((x).^2))/100);

    end
    
    mapindex = Gini;
    index = mean(Gini);

end


if strcmp(method,'KolmogorovSmirnov')
    KS = zeros(dim,1);
    Hs = zeros(dim,1);
    for i = 1:dim
            p = mapdists(:,i);
            [f,x] = ecdf(p);
            [h,p,ks2stat] = kstest2(f,x);
            KS(i) = 1-ks2stat;
            Hs(i) = h;
    end
    mapindex = KS;
    index = mean(Hs);

end

if strcmp(method,'Kurtosis')
    K = zeros(dim,1);
    maxvar = max(var(mapdists,0,2));
    unifKurt = 9/5;
    for i = 1:dim
            p = mapdists(:,i);
            varp = var(p);
            K(i) = (varp./maxvar)*(abs((kurtosis(p)-unifKurt)/max([unifKurt,kurtosis(p)]))+1);
    end

    
    mapindex = K;
    index = mean(K);

end

if strcmp(method,'Rank')
    Rank = zeros(dim,1);
    Hs = zeros(dim,1);
    Ps = zeros(dim,1);
    for i = 1:dim
            p = mapdists(:,i);
            x = rand(length(p),1);
            [ps,h,stats] = ranksum(p,x);
            Rank(i) = stats.ranksum;
            Hs (i) = abs(1-h);
            Ps(i) = 1-ps;

    end

    
    mapindex = Hs;
    index = mean(Hs);

end

if strcmp(method,'kkplot')
    EK = zeros(dim,1);
    for i = 1:dim
            p = mapdists(:,i);
            EK(i) = kurtosis(p)-3;
    end
    unifEK = 6/5*ones(dim,1);
    linmod = fitlm(EK, unifEK)
    plot(linmod)
    mapindex = abs(EK./unifEK);
    index = mean(mapindex);

end
[SortedIndex,idx] = sort(mapindex);

if plotting
    Info.idx = idx;
    Info.SortNormDistance = SortedIndex;
    Map = NicheGeneration(Dimensions, Info, 1, true);
end

met.Metric = index
met.MapMetric = mapindex;
met.SortedNormalizedIndex = SortedIndex;
met.idx = idx;


end
