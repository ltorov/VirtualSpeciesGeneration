function Deformations = HarmonicDeformations(norm_climate_vars, layer_num, deformations_limit, plotting)
% Deformations = HarmonicDeformations(NormalizedClimVar, layer_num, limdef, plotting)
%
% DESCRIPTION: 
%   This function generates harmonic deformations of a given climate variable data set.
%   It randomly selects two principal component axes from the data set and uses them to
%   calculate a radial distance and an angle value for each point in a grid. It then
%   generates harmonic functions for these values and calculates the distance of each
%   point in the grid from the center. Finally, it creates a scatter plot with distance
%   values colored by their corresponding data values.
% 
% INPUTS:
%   norm_climate_vars: normalized climate variable data.
%   layer_num: number of layers.
%   deformations_limit: limit of deformation (default value = 2).
%   plotting: whether to generate a plot (default value = false).
%
% OUTPUTS:
%   Deformations: struct containing distance values for each deformation.

% Set default values for limdef and plotting if they are not provided
if nargin < 3
    deformations_limit = 2;
end
if nargin < 4
    plotting = false;
end

% Calculate the principal component axes of the input data
PCAs = pca(norm_climate_vars');

% Randomly select two principal component axes
ind = randi([1 layer_num], 1, 2);
PCA_axis_a = PCAs(ind(1), :) * norm_climate_vars;
PCA_axis_b = PCAs(ind(2), :) * norm_climate_vars;

% Generate a random boolean value to determine whether to use min or max values for r
boo = randi([0 1], 1);
if boo == 0
    r = (PCA_axis_a - min(PCA_axis_a)) / (max(PCA_axis_a) - min(PCA_axis_a));
else
    r = (PCA_axis_a - max(PCA_axis_a)) / (min(PCA_axis_a) - max(PCA_axis_a));
end

% Calculate the alpha angle value
alpha = ((PCA_axis_b - min(PCA_axis_b)) / (max(PCA_axis_b) - min(PCA_axis_b))) * 2 * pi;

% Generate a grid of sample_num points with random values for r and alpha
sample_num = 1e4;
samples = rand(sample_num, 2) * diag([1, 2 * pi]);

% Generate harmonic functions for each point in the grid and calculate the distance from the center
H = HarmonicFunction(samples, deformations_limit, plotting);
distances = H.distances;

% Create a scattered interpolant object and use it to calculate the distance for each point with the given r and alpha values
F = scatteredInterpolant(samples(:, 1), samples(:, 2), distances);
Z = F(r, alpha);

% If plotting is true, generate a scatter plot with distance values colored by data values and save it as a figure
if plotting
    figure(3)
    scatter(r, alpha, [], Z, 'filled')
    colormap jet
    savefig('ColoredPCA.fig')
end

% Store the distance values in a struct and return it
Deformations.Distances = Z;
end
