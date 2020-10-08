function [t,x,days,K_disc,A,Kfun] = lockdown(data, K0_disc, K0_cont, options)

%
%   [t,x,days,K_disc,A,Kfun] = lockdown(data, K0_disc, K0_cont, options)
%
%   Lockdown è divisa in 3 parti: calcolo i k discreti, fitto i valori
%   discreti per ottenere una funzione continua k(t) inserendo i parametri
%   del fitting polinomiale nell'output A e simulo il modello aggiornato.
%   Inoltre nel caso italiano calcolo e plotto il funzionale L(K).
%
%   INPUTS:
%   data        : struttura che contiene i dati utili come Ibar e Rbar
%   K0_disc     : guess per il calcolo dei K discreti
%   K0_cont     : guess per il calcolo di K continuo
%   options     : struttura contenente i campi che regolano plot e il campo
%                 pnt per aumentare nodi integrazione in minquad_kdiscreti.
%                 Campo ffit per il caso italia per selezionare la
%                 regressione di K(t) tra gaussiana (0) e esponenziale (1).
%
%   OUTPUTS:
%   t           : tempi soluzione del modello simulato
%   x           : soluzione modello simulato
%   days        : giorni (non indici) in cui ho trovato K discreto
%   K_disc      : valori discreti ottenuti del parametro K
%   A           : parametri stimati per fittare i K_disc
%   K_fun       : funzione risultante che fitta i K_disc
%

%global t_u t_c Nass Ibar Rbar beta gamma date K_disc days regione

% recupero i valori che servono
[Nass,Ibar,Rbar] = data.value;
[~,~,t_u,t_c,~] = data.time;
[beta,gamma] = data.parameters;

if isfield(data,'regione')
    ffunz = 0;      % nelle regioni non ho previsto il plot del funzionale
else    % caso italia
    if isfield(options,'ffunz')
        ffunz = options.ffunz;
    else
        ffunz = 1;
    end
end

if nargin == 3
    ffig = 1;
    ssave = 1;
    data(1).ffit = 0;
else
    if isfield(options,'ffig')
        ffig = options.ffig;
    end
    if isfield(options,'ssave')
        ssave = options.ssave;
    end
    if isfield(options,'ffit')
        ffit = options.ffit;
        data(1).ffit = options.ffit;
    end
end

%% Plot Funzionale

% plotto a video il funzionale relativo ai K discreti
% discreti
if ffunz == 1
    plot_funzionale
end

%% 1. stimo i K discreti

% definisco la finestra
window.h   = 1;              % daily time step
window.kl  = 3;
window.kr  = 4;
% non arrivo fino a t_c, altrimenti non avrei una finestra di 7 giorni
% centrata in t_c
kspan      = t_u:1:t_c-window.kr*window.h;     % intervallo ricerca k discreto

pnt = 1;
if isfield(options,'pnt')
    pnt = options.pnt;
end

% per comodità, procedura per i K_disc discreti su kspan su un file a parte
[days, K_disc] = stima_kdiscreti(data,kspan,window,K0_disc,pnt);

data(1).fittingK = days;
data(2).fittingK = K_disc;

%% 2. Fitto i K_disc discreti ottenuti e ricavo k(t)

ffig = 1;
ssave = 0;

[Kfun,A] = stima_kcontinuo(data,K0_cont,ffig,ssave);

%% Simulazione modello oltre il lockdown

% update sistema
SI = @(t,x) [-(beta - x(1)*x(2)/Kfun(t))*x(1)*x(2);
              (beta - x(1)*x(2)/Kfun(t))*x(1)*x(2) - gamma*x(2)];

Jac = @(t,x) [ -beta*x(2) + 2*x(1)*(x(2)^2)/Kfun(t), -beta*x(1) + 2*(x(1)^2)*x(2)/Kfun(t);
                beta*x(2) - 2*x(1)*(x(2)^2)/Kfun(t),  beta*x(1) - 2*(x(1)^2)*x(2)/Kfun(t) - gamma];
options.Jacobian = Jac;

I0 = Ibar(t_u+1); R0 = Rbar(t_u+1); S0 = Nass-I0-R0;
x0 = [S0;I0]/Nass;          % dato iniziale in percentuale

deltatc = 30;
nstep = 500;
if isfield(options,'nstep')
    nstep = options.nstep;
end

[t, x]  = eulerorosenbrock(SI,linspace(t_u,t_c+deltatc,nstep),x0,options);
x(:,3) = ones(length(t),1) - x(:,1) - x(:,2);      % ricavo R per post-processing
x = Nass.*x;                                       % normalizzo

end