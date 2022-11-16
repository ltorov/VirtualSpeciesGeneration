function Deformations = BetaDeformations(NormalizedClimVar,point,NumLayers,n, plotting,boundlo,boundhi)
% Deformations = BetaDeformations(X,point,NumLayers,n, plotting,boundlo,boundhi)
% 
% DESCRIPTION
%   
% 
% REQUIRED INPUTS
%   NormalizedClimVar: A matrix of normalized layers, 
%                    the rows are the number of layers, 
%                    and the columns are the map shaped as a vector. 
%   point: an array of [NumLayers,1] of a chosen initial point.
%   NumLayers: an integer with the number of layers.
%   n: an integer with the length of the layer.
%
% OPTIONAL INPUTS
%   plotting: boolean variable (true, false) to plot the beta functions.
%   boundlo:
%   boundhi:
% 
% OUTPUTS
%   %Deformations: a structure containing
%       -ClimVar:
%       -NewPoint:
%%
    % COMPLETE INPUT
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
    a = boundlo+(boundhi+boundlo).*rand(NumLayers,1);
    b = boundlo+(boundhi+boundlo).*rand(NumLayers,1);
    for i = 1:NumLayers
        if a(i,:)<1
            if b(i,:)>=1
                b(i,:) = rand(1);
            end
        end
        if a(i,:)>=1
            if b(i,:)<1
                b(i,:) = boundlo+(boundhi+boundlo).*rand(1);
            end
        end
    end

    % Warp the space and the initial point using the distributions
    % generated.
    NewPoint = point;
    for i = 1:NumLayers
        ClimVar(i,:) = betacdf(NormalizedClimVar(i,:),a(i),b(i));
        NewPoint(i) = betacdf(point(i),a(i),b(i));
    end

    % Plotting of the beta distributions.
    if plotting == true
        clf
        hold on
        for i =1:NumLayers
            climvar = ClimVar(i,:);
            SortedNormalizedClimVar = sort(NormalizedClimVar,2);
            plot(SortedNormalizedClimVar(i,:),climvar)
        end
    end

    % OUTPUT STORAGE    
    Deformations.ClimVar = ClimVar;
    Deformations.NewPoint = NewPoint;

   
end
