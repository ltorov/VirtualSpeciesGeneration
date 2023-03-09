function ReadInfo = ReadLayers(layerfolder, parallel, nanvalue)
% ReadInfo = ReadLayers(layerfolder, parallel, nanvalue)
%
% DESCRIPTION
% 'ReadLayers' reads multiple environmental layers stored in '.asc' files.
%
% REQUIRED INPUTS 
%   layerfolder: string with the path to the folder containing the layers.
%
% OPTIONAL INPUTS
%   parallel: boolean variable (true, false), perform the file reading in parallel,
%             (default: false).
%   nanvalue: the value to use as NaN (Not a Number) in the layers (default: -9999).
%
% OUTPUTS
%   ReadInfo: a structure containing:
%       -Z: a 3d matrix with the layers information.
%       -R: a matrix of the layers coordinates.
%       -Map: a 2d matrix with the mean value of the layers.
%       -Indicator: a vector with the index of NaN values in the layers.
%       -Dimensions: a vector with the length of the layers 
%                    and the number of layers.  
%       -NormalizedClimVar: A matrix of normalized layers, 
%                    the rows are the number of layers, 
%                    and the columns are the map shaped as a vector.  

% Start timer
tic

% Set default values for optional inputs
if nargin < 2
    parallel = false;
end

if nargin < 3
    nanvalue = -9999;
end

% Set up parallel processing, if requested
if parallel
    parforArg = inf;
else
    parforArg = 0;
end

% Add folder to MATLAB path to read layers
addpath(layerfolder)

% Get list of layer files in folder
LayerDir = dir(layerfolder);
layers = {LayerDir.name};

% Keep only files that start with 'bio'
CompareFileLayers = strncmp('bio', layers, 3);
layers = layers(CompareFileLayers);

% Read the first layer to get its map and coordinates
[Z,R] = readgeoraster(strcat(layerfolder, layers{1}),'CoordinateSystemType','geographic');
N = length(layers);

% Preallocate Z matrix to store all layers
Z = zeros(size(Z, 1), size(Z, 2), N);

disp('----Reading layers----')

% Read each layer and store it in Z matrix
parfor (i = 1 : N, parforArg)
    interZ = readgeoraster(layers{i},'CoordinateSystemType','geographic');
    interZ(interZ <= nanvalue) = NaN;
    Z(:, :, i) = interZ;
end

% Remove folder from MATLAB path
rmpath(layerfolder)

% Compute mean map and remove NaNs from layers
[Rows, Columns, NumLayers] = size(Z);
Map = Z(:, :, 1);
ClimVariables = zeros(NumLayers, Rows * Columns);

parfor (i = 1 : NumLayers, parforArg)
    aux = Z(:, :, i);
    ClimVariables(i, :) = aux(:);
end

nanDetector = sum(ClimVariables);
nanPositions = isnan(nanDetector);

ClimVar = ClimVariables(:,~nanPositions);
Dimension = size(ClimVar,2);
Map(:) = mean(ClimVariables);

% OUTPUT STORAGE
ReadInfo.Indicator = nanPositions;

ReadInfo.Dimensions = [Dimension, NumLayers];
ReadInfo.NormalizedClimVar = normalize(ClimVar, 2, 'range');
ReadInfo.Map = Map;
ReadInfo.Z = Z;
ReadInfo.R = R;

% Stop timer
toc
    
end