function Deformations = HarmonicDeformations(NormalizedClimVar,NumLayers,limdef,plotting)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

if nargin < 3
    limdef = 2;
end
if nargin < 4
    plotting = false;
end


PCAs = pca(NormalizedClimVar');
ind = randi([1 NumLayers],1,2);
PCAa = PCAs(ind(1),:)*NormalizedClimVar; PCAb = PCAs(ind(2),:)*NormalizedClimVar;
boo = randi([0 1],1);
if boo == 0
    r = (PCAa- min(PCAa))/(max(PCAa)-min(PCAa));
else
    r = (PCAa- max(PCAa))/(min(PCAa)-max(PCAa));
end

alpha = ((PCAb- min(PCAb))/(max(PCAb)-min(PCAb)))*2*pi;

%samples = [r ; alpha]';


Nsamples = 1e4;%número de puntos
samples = rand(Nsamples,2) * diag([1,2*pi]);
H = HarmonicFunction(samples, limdef, plotting);
distances = H.distances;


F = scatteredInterpolant(samples(:,1),samples(:,2),distances);
Z = F(r,alpha);

if plotting
    figure(3)
    scatter(r,alpha,[],Z,'filled')
    colormap jet
    savefig('ColoredPCA.fig')
end

Deformations.distances = Z;

end