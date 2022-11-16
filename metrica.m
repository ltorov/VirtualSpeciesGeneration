function met = metrica(maps,method)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
if nargin<2
    method = 'GiniIndex';
end
% method = 'ShannonEntropy' o 'GiniIndex'
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
        Entropy(i) = -sum(p.*log(p+eps));
    end

    
    unifEntropy = -sum(log2(1/k));
    index = Entropy./unifEntropy;

end
if strcmp(method,'GiniIndex')

    Gini = zeros(dim,1);
    for i = 1:dim
        p = mapdists(:,i);
        [f,x] = ecdf(p);
        Gini(i) = (sum(((x-f).^2))/100)/(sum(((x).^2))/100);
        %Gini(i) = (sum(((f).^2))/100)/(sum(((x).^2))/100);
    end

    
    index = Gini;

end


%plot(index,'o')

NormalizedIndex = normalize(index,'range');
[SortedNormalizedIndex,idx] = sort(NormalizedIndex);

index = mean(index);

met.Metric = index
met.NormalizedIndex = NormalizedIndex;
met.SortedNormalizedIndex = SortedNormalizedIndex;
met.idx = idx;


end