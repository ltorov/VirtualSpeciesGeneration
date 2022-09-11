n = 1000;
x1 = linspace(0,1,n)'; x2 = linspace(0,1,n)';
a1 = 5.*rand(n,1); a2 = 5.*rand(n,1);
b1 = 5.*rand(n,1); b2 = 5.*rand(n,1);

for i = 1:n
    if a1(i,:)<1
        if b1(i,:)>=1
            b1(i,:) = rand(1);
        end
    end
    if a1(i,:)>=1
        if b1(i,:)<1
            b1(i,:) = 1 +6.*rand(1);
        end
    end
    if a2(i,:)<1
        if b2(i,:)>=1
            b2(i,:) = rand(1);
        end
    end
    if a2(i,:)>=1
        if b2(i,:)<1
            b2(i,:) = 1 +6.*rand(1);
        end
    end

end
g1 = betacdf(x1,a1,b1); g2 = betacdf(x2,a2,b2);