%{
    Situazione Italia: prelockdown, lockdown
    Alla fine plot dela simulazione del modello in entrambe le fasi.
%}

clear all
close all
clc

% controllo blocchi codice: posso fermare lockdown e riepilogo
lock      = 1;
riepilogo = 1;

%% DATI

% Upload dati protezione civile
tmp = fullfile('..','00 - dpc_data','2020-05-22','dati-andamento-nazionale');
%tmp = fullfile('..','00 - dpc_data','2020-06-30','dati-andamento-nazionale');

[status,result]     = fileattrib(tmp);
path_folder         = result.Name;                  % percorso alla cartella
[date,Ibar,Rbar]    = data_read_dpc(path_folder);
% Attenzione: gli indici in matrice partono da 1 e tm parte da 0
% date(i), Ibar(i), Rbar(i) e' in corrispondenza con t_i-1

% Dati
t_0 = 0;                           % 2020-02-24 iniziale
t_u = 14;                          % 2020-03-09 t finale senza controllo
t_c = length(date)-1;              % ultimo tempo
Nass = 60317000;                   % popolazione italiana istat 11.02.2020
%Nass = 6e9;                       % imbrogliando così arrivo a k~1e-5 (in
                                   % main lockdown

% Creo struttura dati da passare alle function

data(1).value = Nass;
data(2).value = Ibar;
data(3).value = Rbar;

data(1).time = t_0;
data(2).time = t_u;
data(3).time = t_c;
data(4).time = date;

%% Pre-Lockdown [t_0,t_u]: stima beta e gamma prima del Lockdown

% controlo immagini e figure
options.ssens = 0;      % analisi sensitività
options.ffig  = 1;      % stampare figure
options.ssave = 1;      % salvare immagine

options.pnt   = 5;      % piu nodi per migliore risoluzione sistema minquad
K0  = [0.4,0.1];        % guess iniziale ottimizzazione per [beta,gamma]

[tpl, xpl, beta, gamma] = prelock(data, K0, options);

% inserisco i parametri nella struttura "data"
data(1).parameters = beta;
data(2).parameters = gamma;

R_0 = beta/gamma;
table(beta,gamma,R_0)

if lock == 0
    return
else
    clear options   % ripulisco il settaggio figure
end

%% Lockdown

% controlo immagini e figure
options.ffunz = 0;      % stampare funzionale
options.oold  = 0;      % stampo il funzionale nella vecchia maniera
options.ffig  = 1;      % stampare le figura
options.ssave = 1;      % salvare le figure

% 1. stimo i k discreti
K0_disc	= 1e-5;         % guess iniziale
options.pnt	= 100;      % aumento numero nodi integrazione

% 2. Fitto i k discreti ottenuti e ricavo k(t)
a = 1e-6; b = 1e-4; c = 1e-3;       % fitting polinomiale
K0_cont = [a,b,c];                  % guess iniziale

% 3. Simulazione modello oltre il lockdown
options.deltatc = 10;

% simulazione modello durante lockdown e fitting dei k discreti
[tl,xl,A] = lockdown(data, K0_disc, K0_cont, options);

if riepilogo == 0
    return
else
    clear options
end    

%% Figure: riepilogo

ssave = 1;  % salvo la figura

fig = figure();

tt = [tpl;tl]'; ii = [xpl(:,2);xl(:,2)]';

p1 = plot(t_0:t_c,Ibar,'o',...
    'MarkerSize',4,...
    'MarkerEdgeColor','red',...
    'MarkerFaceColor',[1 .6 .6]);

hold on
p2 = plot(tt,ii,'Linewidth',2.5,'color','black');
p2.Color(4) = 0.6;

ax = gca;
ax.XTick = [t_0,t_u,37,67,t_c];
ax.XTickLabel = date([t_0,t_u,37,67,t_c]+1);
ax.XTickLabelRotation = 45;
xline(t_u,':','inizio Lockdown')
%ylim([0 3e5])

box on
legend([p1,p2],'$I_{bar}$','I','Location','Best');
ylabel('casi confermati');
set(gca,'FontSize',12.5)

% imposto latex come inteprete per i grafici
set(groot,...
    'defaulttextinterpreter','latex',...
    'defaultAxesTickLabelInterpreter','latex',...
    'defaultLegendInterpreter','latex');  

if ssave == 1
    exportgraphics(fig,['figure/italia_riepilogo ' num2str(date(t_c+1)) '.pdf'],'ContentType','vector',...
                   'BackgroundColor','none')
end
