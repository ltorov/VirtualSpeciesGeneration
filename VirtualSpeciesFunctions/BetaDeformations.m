function Deformations = BetaDeformations(norm_climate_vars, point, layer_num, layer_len, plotting, lower_bound, upper_bound)
% Deformations = BetaDeformations(normalized_climate_vars, point, layer_num, layer_len, plotting, lower_bound, upper_bound)
% 
% DESCRIPTION:
%   Generates beta distributions and warps a given point according to those
%   distributions.
% 
% REQUIRED INPUTS:
%   norm_climate_vars: A matrix of normalized layers, the rows are
%                            the number of layers, and the columns are the 
%                            map shaped as a vector.
%   point: an array of [layer_num,1] of a chosen initial point.
%   layer_num: an integer with the number of layers.
%   layer_len: an integer with the length of the layer.
%
% OPTIONAL INPUTS:
%   plotting: boolean variable (true, false) to plot the beta functions.
%   lower_bound: the lower bound of the beta distribution.
%   upper_bound: the upper bound of the beta distribution.
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
        lower_bound = 0;
    end
    if nargin < 7
        upper_bound = 5;
    end

    % Generate alpha and beta parameters for NumLayers beta distributions.
    climate_vars = zeros(layer_num, layer_len);
    a = lower_bound + (upper_bound + lower_bound) .* rand(layer_num, 1);
    b = lower_bound + (upper_bound + lower_bound) .* rand(layer_num, 1);
        
    % Ensure a and b values fall within acceptable ranges.
    for i = 1:layer_num
        if a(i, :) < 1 && b(i, :) >= 1
            b(i, :) = rand(1);
        elseif a(i, :) >= 1 && b(i, :) < 1
            b(i, :) = lower_bound + (upper_bound - lower_bound) .* rand(1);
        end
    end

    % Warp the space and the initial point using the distributions
    % generated.
    new_point = point;
    for i = 1:layer_num
        climate_vars(i,:) = betacdf(norm_climate_vars(i,:),a(i),b(i));
        new_point(i) = betacdf(point(i),a(i),b(i));
    end

    % Plotting of the beta distributions if plotting is true.
    if plotting
        clf
        hold on
        for i = 1:layer_num
            climvar = climate_vars(i,:);
            SortedNormalizedClimVar = sort(norm_climate_vars, 2);
            plot(SortedNormalizedClimVar(i,:), climvar)
        end
    end

    % Store the output in a structure.
    Deformations.ClimVar = climate_vars;
    Deformations.NewPoint = new_point;

end