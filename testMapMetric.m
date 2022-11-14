clear;clc;close all
%%
layerfolder='../Tool_landscape/data/Colombia/';%route to folder with ambiental data
Dimensions = ReadLayers(layerfolder);
%%
k = 90;
j = 90;
method = 'harmonic';
tic
for i = j:k
    waitbar((i-j)/(k-j), 'Generating maps')
    Info = InitialPoint(Dimensions, method, false, 'limdef', 3);
    %figure(i);
    Map = NicheGeneration(Dimensions, Info, 1, false);
    MapsH(i) = Map;
end
toc
%%
method = 'coeff';

k = 90;
j = 71;
h = waitbar(0,'Generating maps');
tic
for i = j:k
    waitbar((i-j)/(k-j),h,sprintf('Generating maps'))
    Info = InitialPoint(Dimensions, method, false);
    %figure(i);
    Map = NicheGeneration(Dimensions, Info, 1, false);
    MapsC(i) = Map;
end
close (h)
toc
%%
method = 'beta';

k = 90;
j = 71;
h = waitbar(0,'Generating maps');
tic
for i = j:k
    waitbar((i-j)/(k-j),h,sprintf('Generating maps'))
    Info = InitialPoint(Dimensions, method, false);
    %figure(i);
    Map = NicheGeneration(Dimensions, Info, 1, false);
    MapsB(i) = Map;
end
close (h)
toc
%%

metH = metrica(MapsH);
%%
Info.idx = metH.idx';
Info.SortNormDistance = metH.SortedNormalizedIndex';


MapH = NicheGeneration(Dimensions, Info, 1, true);

%%
metC = metrica(MapsC);


%%
Info.idx = metC.idx';
Info.SortNormDistance = metC.SortedNormalizedIndex';


MapC = NicheGeneration(Dimensions, Info, 1, true);


%%
metB = metrica(MapsB);
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
m=1000;
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
