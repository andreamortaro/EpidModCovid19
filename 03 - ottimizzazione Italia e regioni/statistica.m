% Modelli monte carlo?
function statistica(data)

% parametri modello
alfab = 0;
alfag = 0;

%p = 40;         % percentuale infetti che aggiungo
c = 8.56;
p = 2*(c-1);
M = 30;        % numero simulazioni con dato valore infetti iniziali misurato
B = 10;          % numero simulazioni con dato beta e gamma

data(1).parametersStat = alfab;
data(2).parametersStat = alfag;
data(3).parametersStat = p;
data(4).parametersStat = M;
data(5).parametersStat = B;

%% prelockdown

ffig = 1;
ssave = 0;

options.ffig = ffig;
options.ssave = ssave;

[tPL, ImedioPL,VarmediaPL,I0f,R0f,hist] = prelockStat(data,options);

clear options

%% lockdown

ffig = 1;
ssave = 1;

options.ffig = ffig;
options.ssave = ssave;
options.deltatc = 10;

[tL, ImedioL, VarmediaL,hist] = lockdownStat(data,I0f,R0f,hist,options);
% hist.sim: sulle righe della cell variano i beta e gamma

%% riepilogo

ffig = 1;
ssave = 1;

options.ffig = ffig;
options.ssave = ssave;

tt = [tPL;tL]; ii = [ImedioPL;ImedioL]; vv = [VarmediaPL;VarmediaL];

riepilogoStat(data,tt,ii,vv,hist,options)

return



end
