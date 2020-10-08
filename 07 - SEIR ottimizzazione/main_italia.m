%{
    Situazione Italia: prelockdown, lockdown
    Alla fine plot dela simulazione del modello in entrambe le fasi.
%}

clear all
close all
clc

% controllo blocchi codice: posso fermare lockdown e riepilogo
lock = 1;
riep = 1;

%% DATI

% Upload dati protezione civile
tmp = fullfile('..','00 - dpc_data','2020-06-30','dati-andamento-nazionale');
%tmp = fullfile('..','00 - dpc_data','2020-08-31','dati-andamento-nazionale');

[status,result]     = fileattrib(tmp);
path_folder         = result.Name;                  % percorso alla cartella
[date,Ebar,Ibar,Rbar,totCases]    = data_read_dpc(path_folder);
% Attenzione: gli indici in matrice partono da 1 e tm parte da 0
% date(i), Ibar(i), Rbar(i) e' in corrispondenza con t_i-1

% Dati
t_0 = 0;                           % 2020-02-24 iniziale
t_1 = 6;
t_u = 24;                          % 2020-03-09 t finale senza controllo
%t_u=14;
t_c = length(date)-1;              % ultimo tempo
Nass = 60317000;                   % popolazione italiana istat 11.02.2020
%Nass = 6e9;                       % imbrogliando così arrivo a k~1e-5 (in
                                   % main lockdown

% Creo struttura dati da passare alle function

data(1).value = Nass;
data(2).value = Ibar;
data(3).value = Rbar;
data(4).value = totCases;
data(5).value = Ebar;

data(1).time = t_0;
data(2).time = t_1;
data(3).time = t_u;
data(4).time = t_c;
data(5).time = date;

%% Pre-Lockdown [t_0,t_u]: stima beta e gamma prima del Lockdown

% controlo immagini e figure
options.ffig  = 1;      % stampare figure
options.ssave = 1;      % salvare immagine

options.pnt   = 5;      % piu nodi per migliore risoluzione sistema minquad
K0  = [0.2,0.05,0.1];        % guess iniziale ottimizzazione per [beta,gamma]

[tpl, xpl, beta, gamma,mu] = prelock(data, K0, options);

% inserisco i parametri nella struttura "data"
data(1).parameters = beta;
data(2).parameters = gamma;
data(3).parameters = mu;

R_0 = beta/gamma;
table(beta,gamma,mu,R_0)

if lock == 0
    return
else
    clear options   % ripulisco il settaggio figure
end

%% Lockdown

% controlo immagini e figure
options.ffig  = 1;      % stampare le figura
options.ssave = 1;      % salvare le figure

% 1. stimo i k discreti
K0_disc	= 1e-3;         % guess iniziale
options.pnt	= 1000;      % aumento numero nodi integrazione

% 2. Fitto i k discreti ottenuti e ricavo k(t)
options.ffit = 0;       % 0 gaussiana, 1 esponenziale

switch options.ffit
    case 0
        a = 0.01; b = 60; c = 30;     % guassiana
    case 1
        a = 1; b = 0.05; c = 0.01;       % exp^3
end

K0_cont = [a,b,c];                  % guess iniziale

% 3. Simulazione modello oltre il lockdown
options.nstep = 2000;

% simulazione modello durante lockdown e fitting dei k discreti
[tl,xl,days,K_disc,A,Kfun] = lockdown(data, K0_disc, K0_cont, options);
data(1).Kvalue = days;
data(2).Kvalue = K_disc;
data(3).Kvalue = Kfun;
data(4).Kvalue = A;

if riep == 0
    return
else
    clear options
end    

%% Figure: riepilogo

options.ssave = 1;

tt = [tpl;tl]'; ii = [xpl(:,3);xl(:,3)]'; ee = [xpl(:,2);xl(:,2)]';

riepilogo(data,tt,ii,ee,options);

