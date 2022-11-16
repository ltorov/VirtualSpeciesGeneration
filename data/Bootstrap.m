%Estimar la media
%y = exprnd(5,200,1); %exponencial
y = 10+2*randn(200,1);
y = [y;20]; %outlier
m = bootstrp(1000,@mean,y);

esperanzaBootstrap = mean(m)
desviacionBootstrap = std(m)
mediaMuestral = mean(y)
%hist(m)
%h:hipotesis, p: pvalue, ci:confidence interval 95%
[h p ci] = ttest(y)
%bootstrap confidence interval
CIB = [prctile(m,2.5) prctile(m,97.5)] 
%%
%Estimar la mediana
%Recordar que la mediana es sesgada, entonces no va a dar igual
y = exprnd(5,100,1); %exponencial
m = bootstrp(1000,@median,y);
esperanzaBootstrap = mean(m)
desviacionBootstrap = std(m)
medianaMuestral = median(y)
%hist(m)
%bootstrap confidence interval
CIB = [prctile(m,2.5) prctile(m,97.5)] 
%%
%Estimar la varianza
%Recordar que la mediana es sesgada, entonces no va a dar igual
y = exprnd(5,100,1); %exponencial
m = bootstrp(1000,@var,y);
esperanzaBootstrap = mean(m)
desviacionBootstrap = std(m)
medianaMuestral = var(y)
hist(m)
%bootstrap confidence interval
CIB = [prctile(m,2.5) prctile(m,97.5)] 