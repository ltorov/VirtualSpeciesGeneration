clear;clc;close all
%%
layerfolder='../VirtualSpeciesGeneration/data/layers/'
Dimensions = ReadLayers(layerfolder);
%%
k = 100;
j = 1;
method = 'harmonic';
h = waitbar(0,'Generating maps');
tic
for i = j:k
    waitbar((i-j)/(k-j),h, ['Generating map ' num2str(i)])
    Info = InitialPoint(Dimensions, method, false, 'limdef', 5);
    %figure(i);
    Map = NicheGeneration(Dimensions, Info, 1, false);
    MapsH(i) = Map;
end
close(h)
toc
%%
method = 'coeff';

k = 100;
j = 1;
h = waitbar(0,'Generating maps');
tic
for i = j:k
    waitbar((i-j)/(k-j),h, ['Generating map ' num2str(i)])
    Info = InitialPoint(Dimensions, method, false);
    %figure(i);
    Map = NicheGeneration(Dimensions, Info, 1, false);
    MapsC(i) = Map;
end
close (h)
toc
%%
method = 'beta';

k = 100;
j = 1;
h = waitbar(0,'Generating maps');
tic
for i = j:k
    waitbar((i-j)/(k-j),h, ['Generating map ' num2str(i)])
    Info = InitialPoint(Dimensions, method, false);
    %figure(i);
    Map = NicheGeneration(Dimensions, Info, 1, false);
    MapsB(i) = Map;
end
close (h)
toc
%%
%method = 'LorenzCurve','KolmogorovSmirnov','ShannonEntropy','Rank','kkplot'
method = 'LorenzCurve';
plotting = false;
metH = metrica(MapsH, Dimensions, method, plotting);
%%
Info.idx = metH.idx';
Info.SortNormDistance = metH.SortedNormalizedIndex';


MapH = NicheGeneration(Dimensions, Info, 1, true);

%%
%method = 'LorenzCurve','KolmogorovSmirnov','ShannonEntropy','Rank','kkplot'
method = 'KolmogorovSmirnov';
plotting = false;
metC = metrica(MapsC, Dimensions, method, plotting);


%%
Info.idx = metC.idx';
Info.SortNormDistance = metC.SortedNormalizedIndex';


MapC = NicheGeneration(Dimensions, Info, 1, true);


%%
%method = 'LorenzCurve','KolmogorovSmirnov','ShannonEntropy','Rank','kkplot'
method = 'Rank';
plotting = false;
metB = metrica(MapsB, Dimensions, method, plotting);


%%
Info.idx = metB.idx';
Info.SortNormDistance = metB.SortedNormalizedIndex';


MapB = NicheGeneration(Dimensions, Info, 1, true);


%%
m=1000000;
x = zeros(m,1);
k = 10;
for i =1:m
    p = rand(k,1);
    x(i) = -sum(p.*log(p+eps))/k;
end

mean(x)

%%
m = 1000;
x = zeros(m,1);
k = 1000;
for i =1:m
    p = rand(k,1);
    x(i) = -sum(p.*log(p+eps))/k;
end

mean(x)

%%
hist(p)
hold on
ksdensity(p)


%%
maps = MapsB;
k = length(maps);


dim = length(maps(1).NormDistance);

mapdists = zeros(k,dim);

for i = 1:k
    mapdists(i,:) = maps(i).NormDistance;
end


KS = zeros(dim,1);
Hs = zeros(dim,1);
for i = 1:dim
        p = mapdists(:,i);
        [f,x] = ecdf(p);
        [h,p,ks2stat] = kstest2(f,x);
        KS(i) = 1-ks2stat;
        Hs(i) = h;
end

    
met = mean(Hs)

%%

[SortedNormalizedIndex,idx] = sort(KS);
Info.idx = idx;
Info.SortNormDistance = SortedNormalizedIndex';


MapC = NicheGeneration(Dimensions, Info, 1, true);
%%
