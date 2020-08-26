%{

    Situazione Italia: prelockdown, lockdown

    Alla fine plot dela simulazione del modello in entrambe le fasi.

%}

close all
clear
clc

global t_0 t_u t_c Nass Ibar Rbar date beta gamma

% controllo blocchi codice: posso attivare lockdown e riepilogo
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

% DATI:
t_0 = 0;                           % 2020-02-24
t_u = 14;                          % 2020-03-09 t finale senza controllo
t_c = length(date)-1;              % decremento perche' parto da 0
Nass = 60317000;                   % popolazione italiana istat 11.02.2020
%Nass = 6e9;                       % imbrogliando così arrivo a k~1e-5 (in
                                   % main lockdown
                                   
%% Pre-Lockdown

% controlo immagini e figure
ssens = 0;      % analisi sensitività
ffig  = 1;      % stampare figure
ssave = 0;      % salvare immagine

K0  = [0.5,0.1];                    % guess iniziale per [beta,gamma]
pnt = 5;                            % piu nodi per migliore risoluzione sistema minquad

[tpl, xpl, beta, gamma] = main_prelock(K0, pnt, ssens, ffig, ssave);


if lock == 0
    return
end    

%% Lockdown

ffunz = 0;  % stampare funzionale
oold  = 0;  % stampo il funzionale nella vecchia maniera
ffig  = 1;  % stampare le figura
ssave = 1;  % salvare le figure

% 1. stimo i K discreti

k0_c	= 1e-5;          % guess iniziale (ottengo sempre gli stessi k_c)
pnt     = 100;           % aumento numero nodi integrazione

% 2. Fitto i k_c discreti ottenuti e ricavo k(t)

a = 1e-6; b = 1e-4; c = 1e-3;       % fitting polinomiale
kguess = [a,b,c];                   % guess iniziale

% 3. Simulazione modello oltre il lockdown
deltatc = 10;

% simulazione con ER del modello durante il lockdwon con i parametri stimati
% A = [a,b,c] parametri per il fitting polinomiale
[tl,xl,A] = main_lockdown(k0_c, pnt, kguess, deltatc, ffunz, oold, ffig, ssave);

if riepilogo == 0
    return
end    

%% Figure: riepilogo

ssave = 1;  % salvo la figura

fig = figure();

tt = [tpl;tl]';
ii = [xpl(:,2);xl(:,2)]';

%plot(tt,ii,'SeriesIndex',3,'Linewidth',2)
plot(tt,ii,'Linewidth',2,'color',[0 0 0]+0.3)

hold on
plot(t_0:t_c,Ibar,'o',...
    'MarkerSize',3,...
    'MarkerEdgeColor','red',...
    'MarkerFaceColor',[1 .6 .6])

ax = gca;
ax.XTick = [t_0,t_u,37,67,t_c];
ax.XTickLabel = date([t_0,t_u,37,67,t_c]+1);
ax.XTickLabelRotation = 45;
xline(t_u,':','inizio Lockdown')
%ylim([0 3e5])      % come albi

box on
legend('I','$I_{bar}$','Location','Best');
ylabel('casi confermati');
set(gca,'FontSize',12.5)

% imposto latex come inteprete per i grafici
set(groot,...
    'defaulttextinterpreter','latex',...
    'defaultAxesTickLabelInterpreter','latex',...
    'defaultLegendInterpreter','latex');  

if ssave == 1
    exportgraphics(fig,['italia_riepilogo ' num2str(date(t_c+1)) '.pdf'],'ContentType','vector',...
                   'BackgroundColor','none')
end
