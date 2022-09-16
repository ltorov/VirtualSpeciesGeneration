function deformations = BetaDeformations(X,NumLayers,n, plotting)
    if nargin < 4
        plotting = false;
    end
    c = 5;
    G = zeros(NumLayers, n);
    a = c.*rand(NumLayers,1);
    b = c.*rand(NumLayers,1);
    for i = 1:NumLayers
        if a(i,:)<1
            if b(i,:)>=1
                b(i,:) = rand(1);
            end
        end
        if a(i,:)>=1
            if b(i,:)<1
                b(i,:) = 1 +(c+1).*rand(1);
            end
        end
    end
    for i = 1:NumLayers
        G(i,:) = betacdf(X(i,:),a(i),b(i));
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
    deformations = G;
end