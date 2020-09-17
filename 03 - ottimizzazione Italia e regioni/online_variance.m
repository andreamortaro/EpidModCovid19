% Algoritmo Welford per la varianza=sigma^2 (varianza interattiva)
% scelto per evitare cancellazione catastrofica
function variance = online_variance(data)
n     = 0;      % numero dati
media = 0;
M2    = 0;
for x = data
     n = n + 1;                         % aggiungo un dato
     delta = x - media;
     media = media + delta/n;           % costruisco la media empirica ad ogni passo
     M2 = M2 + delta*(x - media);
end
variance = M2/(n - 1);
end