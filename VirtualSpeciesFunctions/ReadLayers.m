function ReadInfo = ReadLayers(layerfolder, parallel,nanvalue)
% ReadInfo = ReadLayers(layerfolder, parallel)
%
% DESCRIPTION
% 'ReadLayers' read multiple environmental layers stored in '.asc' files
%
% REQUIRED INPUTS 
%   layerfolder: string with the layers path 
%
% OPTIONAL INPUTS
%   parallel: bolean variable (true, false), perform the file reading in parallel,
%             (default: false)
%
% OUTPUTS
%   ReadInfo: a structure containing:
%       -Z: a 3d matrix with the layers information
%       -R: a matrix of the layers coordinates
%       -Map: a 2d matrix with the mean value of the layers
%       -Indicator: a vector with the index of NaN values in the layers
%       -Dimensions: a vector with the length of the layers 
%                    and the number of layers  
%       -NormalizedClimVar: A matrix of normalized layers, 
%                    the rows are the number of layers, 
%                    and the columns are the map shaped as a vector.  
%%
    tic
    
    if nargin < 2
        parallel = false;
    end

     if nargin < 3
        nanvalue = -9999;
    end

    if parallel
        parforArg = inf;
    else
        parforArg = 0;
    end
    
    %Layer lecture
    addpath(layerfolder)
    
    LayerDir = dir(layerfolder);
    layers = {LayerDir.name};
    CompareFileLayers = strncmp('bio', layers, 3);
    layers = layers(CompareFileLayers);
    
    %Read map and it coordinates
    [Z,R] = readgeoraster(strcat(layerfolder, layers{1}),'CoordinateSystemType','geographic');
    N = length(layers);
    Z = zeros(size(Z, 1), size(Z, 2), N);

    disp('----Reading layers----')
    
    parfor (i = 1 : N, parforArg)
        interZ = readgeoraster(layers{i},'CoordinateSystemType','geographic');
        interZ(interZ <= nanvalue) = NaN;
        Z(:, :, i) = interZ;
    end
    
    %Create a new data set 
    [Rows, Columns, NumLayers] = size(Z);
    Map = Z(:, :, 1);
    ClimVariables = zeros(NumLayers, Rows * Columns);
    
    parfor (i = 1 : NumLayers, parforArg)
        aux = Z(:, :, i);
        ClimVariables(i, :) = aux(:);
    end
    
    rmpath(layerfolder)
    
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
    
    toc
    
end