function Deformations = BetaDeformations(NormalizedClimVar,point,NumLayers,n, plotting)
% [G, NewPoint] = BetaDeformations(X,point,NumLayers,n, plotting)
% 
% DESCRIPTION
%   
% 
% REQUIRED INPUTS
%   NormalizedClimVar:
%   point:
%   NumLayers:
%   n:
%   plotting:
% 
% OUTPUTS
%   %Deformations: a structure containing
%       -ClimVar:
%       -NewPoint:
%%
    if nargin < 5
        plotting = false;
    end
    boundlo = 0;
    boundhi = 5;
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
    NewPoint = point;
    for i = 1:NumLayers
        ClimVar(i,:) = betacdf(NormalizedClimVar(i,:),a(i),b(i));
        NewPoint(i) = betacdf(point(i),a(i),b(i));
    end
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