function Deformations = HarmonicDeformations(NormalizedClimVar,NumLayers,limdef,plotting)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

if nargin < 3
    limdef =2;
end
if nargin <4
    plotting = false;
end


PCAs = pca(NormalizedClimVar');
ind = randi([0 NumLayers],1,2);
PCAa = PCAs(ind(1),:)*NormalizedClimVar; PCAb = PCAs(ind(2),:)*NormalizedClimVar;
boo = randi([0 1],1);
if boo == 0
    r = (PCAa- min(PCAa))/(max(PCAa)-min(PCAa));
else
    r = (PCAa- max(PCAa))/(min(PCAa)-max(PCAa));
end

alpha = ((PCAb- min(PCAb))/(max(PCAb)-min(PCAb)))*2*pi;

samples = [r ; alpha]';

H = HarmonicFunction(samples, limdef, plotting);

Deformations.distances = H.distances;
%plot(r,alpha,'o')
end