function Deformations = BetaDeformations(normalized_climate_vars, point, layer_num, layer_len, plotting, boundlo, boundhi)
% Deformations = BetaDeformations(normalized_climate_vars, point, layer_num, layer_len, plotting, boundlo, boundhi)
% 
% DESCRIPTION:
%   Generates beta distributions and warps a given point according to those
%   distributions.
% 
% REQUIRED INPUTS:
%   normalized_climate_vars: A matrix of normalized layers, the rows are
%                            the number of layers, and the columns are the 
%                            map shaped as a vector.
%   point: an array of [layer_num,1] of a chosen initial point.
%   layer_num: an integer with the number of layers.
%   layer_len: an integer with the length of the layer.
%
% OPTIONAL INPUTS:
%   plotting: boolean variable (true, false) to plot the beta functions.
%   boundlo: the lower bound of the beta distribution.
%   boundhi: the upper bound of the beta distribution.
% 
% OUTPUTS:
%   Deformations: a structure containing:
%       - ClimVar: a matrix of the beta values for each layer.
%       - NewPoint: the warped point.

    % Set default values for optional inputs.
    if nargin < 5
        plotting = false;
    end
    if nargin < 6
        boundlo = 0;
    end
    if nargin < 7
        boundhi = 5;
    end

    % Generate alpha and beta parameters for NumLayers beta distributions.
    ClimVar = zeros(layer_num, layer_len);
    a = boundlo + (boundhi + boundlo) .* rand(layer_num, 1);
    b = boundlo + (boundhi + boundlo) .* rand(layer_num, 1);
    
    % Ensure a and b values fall within acceptable ranges.
    for i = 1:layer_num
        if a(i, :) < 1
            if b(i, :) >= 1
                b(i, :) = rand(1);
            end
        end
        if a(i, :) >= 1
            if b(i, :) < 1
                b(i, :) = boundlo + (boundhi + boundlo) .* rand(1);
            end
        end
    end

    % Warp the space and the initial point using the distributions
    % generated.
    NewPoint = point;
    for i = 1:layer_num
        ClimVar(i,:) = betacdf(normalized_climate_vars(i,:),a(i),b(i));
        NewPoint(i) = betacdf(point(i),a(i),b(i));
    end

    % Plotting of the beta distributions if plotting is true.
    if plotting == true
        clf
        hold on
        for i =1:layer_num
            climvar = ClimVar(i,:);
            SortedNormalizedClimVar = sort(normalized_climate_vars,2);
            plot(SortedNormalizedClimVar(i,:),climvar)
        end
    end

    % Store the output in a structure.
    Deformations.ClimVar = ClimVar;
    Deformations.NewPoint = NewPoint;

end