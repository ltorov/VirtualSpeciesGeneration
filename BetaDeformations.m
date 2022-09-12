NumLayers = 19;
n = 1000;
x1 = linspace(0,1,n)'; x2 = linspace(0,1,n)';
X = zeros(NumLayers, n);

for i = 1:NumLayers
    X(i,:) = linspace(0,1,n);
end

G = zeros(NumLayers, n);

a = 5.*rand(NumLayers,1);
b = 5.*rand(NumLayers,1);

for i = 1:NumLayers
    if a(i,:)<1
        if b(i,:)>=1
            b(i,:) = rand(1);
        end
    end
    if a(i,:)>=1
        if b(i,:)<1
            b(i,:) = 1 +6.*rand(1);
        end
    end
end

for i = 1:NumLayers
    G(i,:) = betacdf(X(i,:),a(i),b(i));
end

%%
hold on
for i =1:NumLayers
    g = G(i,:);
    plot(X(i,:),g)
end
