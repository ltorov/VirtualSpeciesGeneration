function Deformations = HarmonicDeformations(NormalizedClimVar,NumLayers)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

PCAs = pca(NormalizedClimVar');
ind = randi([0 NumLayers],1,2);
PCAa = PCAs(ind(1),:)*NormalizedClimVar; PCAb = PCAs(ind(2),:)*NormalizedClimVar;


ind = randi([0 1],1);

if ind == 0
    r = (PCAa- min(PCAa))/(max(PCAa)-min(PCAa));
else
    r = (PCAa- max(PCAa))/(min(PCAa)-max(PCAa));
end

alpha = ((PCAb- min(PCAb))/(max(PCAb)-min(PCAb)))*2*pi;

Deformations = [r ; alpha];
end