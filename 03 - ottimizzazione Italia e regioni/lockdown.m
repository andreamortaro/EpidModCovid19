function [t,x,days,K_disc,A,Kfun] = lockdown(data, K0_disc, K0_cont, options)

%
%   [t,x,A] = lockdown(data, K0_disc, kguess, otpions)
%
%   Lockdown è divisa in 3 parti: calcolo i k discreti, fitto i valori
%   discreti per ottenere una funzione continua k(t) inserendo i parametri
%   del fitting polinomiale nell'output A e simulo il modello aggiornato.
%
%   INPUTS:
%   data        : struttura che contiene i dati utili come Ibar e Rbar
%   K0_disc     : guess per il calcolo dei k discreti
%   options     : struttura contenente i campi che regolano plot e il campo
%                 pnt per aumentare nodi integrazione in minquad_kdiscreti
%
%   OUTPUTS:
%   t           : tempi soluzione del modello simulato
%   x           : soluzione modello simulato
%   A           : parametri stimati per fittare i K_disc
%   days        : giorni (non indici) in cui ho trovato k discreto
%   K_disc      : il valore di k discreto
%

%global t_u t_c Nass Ibar Rbar beta gamma date K_disc days regione

% recupero i valori che servono
[Nass,Ibar,Rbar] = data.value;
[~,~,t_u,t_c,~] = data.time;
[beta,gamma] = data.parameters;

if isfield(data,'regione')
    ffunz = 0;      % nelle regioni non ho previsto il plot del funzionale
else    % in questo caso son sicuro di essere in italia
    if isfield(options,'ffunz')
        ffunz = options.ffunz;
    else
        ffunz = 1;
    end
end

if nargin == 3
    ffig = 1;
    ssave = 1;
else
    if isfield(options,'ffig')
        ffig = options.ffig;
    end
    if isfield(options,'ssave')
        ssave = options.ssave;
    end
end

%% Plot Funzionale

if ffunz == 1
    plot_funzionale
end

%% 1. stimo i K discreti

% definisco la finestra
window.h   = 1;              % daily time step
window.kl  = 3;
window.kr  = 4;

% non arrivo a t_c altrimenti in t_c non ho una finestra di 7 giorni
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

deltatc = 0;
nstep = 500;
if isfield(options,'deltatc')
    deltatc = options.deltatc;
    nstep = options.nstep;
end

[t, x]  = eulerorosenbrock(SI,linspace(t_u,t_c+deltatc,nstep),x0,options);
x(:,3) = ones(length(t),1) - x(:,1) - x(:,2);      % ricavo R per post-processing
x = Nass.*x;                                       % normalizzo

end