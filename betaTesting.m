NumLayers = 1;
n = 306593;
X = zeros(NumLayers, n);
for i = 1:NumLayers
    X(i,:) = linspace(0,1,n);
end
BetaDeformations(X,NumLayers,n,true);


%%
NumLayers = Dimensions.Dimensions(2);
NormalizedClimVar = (Dimensions.NormalizedClimVar);
NormalizedClimVar = sort(NormalizedClimVar,2);
n = size(NormalizedClimVar);
n = n(2);
NormalizedClimVar  = BetaDeformations(NormalizedClimVar,NumLayers,n,true);
size(NormalizedClimVar)
%%
point = rand(NumLayers, 1);

point = BetaDeformations(point,NumLayers,1,false);
size(point);
%%
Distance = zeros(1, n);

for i = 2 : n-1
        Distance(i) = norm(point - NormalizedClimVar(:, i))...
                      * (2 - corr2(point, NormalizedClimVar(:, i)));
end
%%

NormDistance = 1 - normalize(Distance, 2, 'range');

%%
[SortNormDistance, idx] = sort(NormDistance, 2, 'descend');

%%
for i =2:5
    i
end
