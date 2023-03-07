function Deformations = BetaDeformations(NormalizedClimVar, point, NumLayers, n, plotting, boundlo, boundhi)
% BETADEFORMATIONS calculates the beta deformations of a given point
%
% INPUTS:
%   NormalizedClimVar: Matrix of normalized layers, where rows represent 
%                      the number of layers and columns represent the 
%                      map shaped as a vector.
%   point: An array of [NumLayers,1] representing a chosen initial point.
%   NumLayers: Integer specifying the number of layers.
%   n: Integer specifying the length of the layer.
%
% OPTIONAL INPUTS:
%   plotting: Boolean variable (true, false) to plot the beta functions.
%   boundlo: Lower bound for alpha and beta parameters. Default is 0.
%   boundhi: Upper bound for alpha and beta parameters. Default is 5.
%
% OUTPUT:
%   Deformations: Structure containing
%       - ClimVar: Beta distributions generated using alpha and beta parameters.
%       - NewPoint: Point after beta deformations.

    % Check if optional inputs are specified
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
    ClimVar = zeros(NumLayers, n);
    a = boundlo + (boundhi + boundlo) .* rand(NumLayers, 1);
    b = boundlo + (boundhi + boundlo) .* rand(NumLayers, 1);
    for i = 1:NumLayers
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

    % Warp the space and the initial point using the distributions generated.
    NewPoint = point;
    for i = 1:NumLayers
        ClimVar(i, :) = betacdf(NormalizedClimVar(i, :), a(i), b(i));
        NewPoint(i) = betacdf(point(i), a(i), b(i));
    end

    % Plot the beta distributions if plotting is enabled.
    if plotting == true
        clf
        hold on
        for i =1:NumLayers
            climvar = ClimVar(i, :);
            SortedNormalizedClimVar = sort(NormalizedClimVar, 2);
            plot(SortedNormalizedClimVar(i, :), climvar)
        end
    end

    % Store the output
    Deformations.ClimVar = ClimVar;
    Deformations.NewPoint = NewPoint;
end
