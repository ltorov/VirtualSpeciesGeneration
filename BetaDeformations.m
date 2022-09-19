function [G, NewPoint] = BetaDeformations(X,point,NumLayers,n, plotting)
    if nargin < 5
        plotting = false;
    end
    c = 5;
    boundlo = 0;
    boundhi = 5;
    G = zeros(NumLayers, n);
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
        G(i,:) = betacdf(X(i,:),a(i),b(i));
        NewPoint(i) = betacdf(point(i),a(i),b(i));
    end
    if plotting == true
        clf
        hold on
        for i =1:NumLayers
            g = G(i,:);
            Xsorted = sort(X,2);
            plot(Xsorted(i,:),g)
        end
    end
end