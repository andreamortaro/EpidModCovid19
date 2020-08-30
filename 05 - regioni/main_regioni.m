%{
    Situazione regioni: prelockdown, lockdown
    Alla fine plot dela simulazione del modello in entrambe le fasi.
%}

clear all
close all
clc

% controllo blocchi codice: posso fermare lockdown e riepilogo
lock      = 1;
riepilogo = 1;

% Upload dati protezione civile
tmp = fullfile('..','00 - dpc_data','2020-05-22','dati-regioni');
[status,result] = fileattrib(tmp);
path_folder = result.Name;              % percorso alla cartella
[reg_label,date,Ibar_allreg,Rbar_allreg] = data_read_dpc_regioni(path_folder);
% Attenzione: poiche' gli indici in matrice partono da 1 e tm parte da 0
% date(i), Ibar(i), Rbar(i) e' in corrispondenza con t_i-1

N_reg = (1:21)';
T= table(N_reg,reg_label) %#ok<NOPTS>

% le regioni con pochi contagli ti portano a K = K0.
i_reg = input("Scegliere n per la corrispondente regione: ");

if(i_reg==3 || i_reg==18)
   warning('foo:bar','N non é corretto per Trento e Bolzano.\n...Press any key to continue')
   pause()
end

% Recupero Nass, Ibar e Rbar della regione i_reg selezionata
tmp = 'tavola_bilancio_mensile_2019_fixed.csv';
path_istat  = fullfile(pwd,tmp);
Nass_reg    = readmatrix(path_istat,'Range',[5 10 25 10],'ExpectedNumVariables',21);

regione = reg_label(i_reg);
Nass    = Nass_reg(i_reg);                     % popolazione regione i_reg
Ibar    = Ibar_allreg{i_reg,:};
Rbar    = Rbar_allreg{i_reg,:};
% NOTA: per trento e bolzano ho inserito un unico valore Nass

% DATI:
t_0 = 0;                    % 2020-02-24
t_u = 14;                   % 2020-03-09 t finale senza controllo
t_c = length(date)-1;       % decremento perche' parto da 0

% Creo struttura dati da passare alle function

data(1).value = Nass;
data(2).value = Ibar;
data(3).value = Rbar;

data(1).time = t_0;
data(2).time = t_u;
data(3).time = t_c;
data(4).time = date;

data(1).regione = regione;

%% Pre-Lockdown [t_0,t_u]: stima beta e gamma prima del Lockdown

% controlo immagini e figure
options.ffig  = 1;      % stampare figure
options.ssave = 1;      % salvare immagine

options.pnt   = 100;     % piu nodi per migliore risoluzione sistema minquad
K0 = [0.5,0.05];        % guess iniziale [beta,gamma]

[tpl,xpl,beta,gamma] = regioni_prelock(data,K0,options);

% inserisco i parametri nella struttura "data"
data(1).parameters = beta;
data(2).parameters = gamma;

R_0 = beta/gamma;
table(regione,beta,gamma,R_0)

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
K0_disc = 1;               % guess iniziale (ottengo sempre gli stessi k_c)
options.pnt = 10;          % aumento numero nodi

% 2. Fitto i k discreti ottenuti e ricavo k(t)

switch regione
    case 'Veneto'
        a = 1e-2; b = 45; c = 37;           % fitting gaussiana
    case {'Emilia-Romagna','Piemonte'}
        a = 1e-3; b = 50; c = 30;
    otherwise
        a = 1e-6; b = 1e-4; c = 1e-3;       % fitting polinomiale
end
K0_cont = [a,b,c];                  % guess iniziale
% cambia in regioni_lockdown e minquad_kcontinuo

% 3. Simulazione modello oltre il lockdown
options.nstep = 1000;
options.deltatc = 10;

% trovo i k_c discreti nell'intervallo [t_u, t_c]
% mi salvo i K_disc per capire il fitting da fare
[tl,xl,A,days,K_disc] = regioni_lockdown(data,K0_disc,K0_cont,options);

if riepilogo == 0
    return
else
    clear options
end    

%% Figure: riepilogo

ssave = 1;  % salvo la figura

fig = figure();

%tt = [tpl;tl]'; ii = [xpl(:,2);xl(:,2)]';

p1 = plot(t_0:t_c,Ibar,'o',...
    'MarkerSize',4,...
    'MarkerEdgeColor','red',...
    'MarkerFaceColor',[1 .6 .6]);

hold on
p2 = plot([tpl;tl]',[xpl(:,2);xl(:,2)]','Linewidth',2.5,'color','black');
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
title(char(regione),'FontSize',13);

set(gca,'FontSize',12.5)

% imposto latex come inteprete per i grafici
set(groot,...
    'defaulttextinterpreter','latex',...
    'defaultAxesTickLabelInterpreter','latex',...
    'defaultLegendInterpreter','latex');  

if ssave == 1
        exportgraphics(fig,...
                       ['figure/' char(regione) '_riepilogo' num2str(date(t_c+1)) '.pdf'],...
                        'ContentType','vector',...
                        'BackgroundColor','none')
end
