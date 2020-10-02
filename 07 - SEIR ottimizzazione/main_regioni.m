%{
    Situazione regioni: prelockdown, lockdown
    Alla fine plot dela simulazione del modello in entrambe le fasi.
%}

clear all
close all
clc

% controllo blocchi codice: posso fermare lockdown e riepilogo
lock = 1;
riep = 1;
stat = 0;

% Upload dati protezione civile
tmp = fullfile('..','00 - dpc_data','2020-06-30','dati-regioni');

[status,result] = fileattrib(tmp);
path_folder = result.Name;              % percorso alla cartella
[reg_label,date,Ibar_allreg,Rbar_allreg,totCases_allreg] = data_read_dpc_regioni(path_folder);
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
totCases= totCases_allreg{i_reg,:};
% NOTA: per trento e bolzano ho inserito un unico valore Nass

% DATI:
t_0 = 0;                    % 2020-02-24
t_1 = 6;
t_u = 14;                   % 2020-03-09 t finale senza controllo
t_c = length(date)-1;       % decremento perche' parto da 0

% Creo struttura dati da passare alle function

data(1).value = Nass;
data(2).value = Ibar;
data(3).value = Rbar;
data(4).value = totCases;

data(1).time = t_0;
data(2).time = t_1;
data(3).time = t_u;
data(4).time = t_c;
data(5).time = date;

data(1).regione = regione;

%% Pre-Lockdown [t_0,t_u]: stima beta e gamma prima del Lockdown

% controlo immagini e figure
options.ffig  = 1;      % stampare figure
options.ssave = 1;      % salvare immagine
options.prelockopt = 0; % minimizzazione in Feb 24 fino a Mar 9


options.pnt   = 100;    % piu nodi per migliore risoluzione sistema minquad
%K0 = [0.3,0.1];           % guess iniziale [beta,gamma]
K0  = [0.3088,0.0495];

[tpl,xpl,beta,gamma] = prelock(data,K0,options);

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
        %a = 0.3; b = 0.04; c = 0.02;         % fitting gaussiana
          a = .004; b = 50; c = 30;
    case 'Lombardia'
         %a = 0.7; b = .02; c = .02;        % exp
         a = .08; b = 50; c = 50;
    case 'Emilia-Romagna'
          a = .004; b = 50; c = 30;
    otherwise
        a = 1e-6; b = 1e-4; c = 1e-3;       % fitting polinomiale
end
K0_cont = [a,b,c];                  % guess iniziale

% 3. Simulazione modello oltre il lockdown
options.nstep = 1000;
options.deltatc = 10;

% trovo i k_c discreti nell'intervallo [t_u, t_c]
% mi salvo i K_disc per capire il fitting da fare
%[tl,xl,A,days,K_disc] = regioni_lockdown(data,K0_disc,K0_cont,options);
[tl,xl,days,K_disc,A,Kfun] = lockdown(data, K0_disc, K0_cont, options);
data(1).Kvalue = days;
data(2).Kvalue = K_disc;
data(3).Kvalue = Kfun;

if riep == 0
    return
else
    clear options
end    
%% Figure: riepilogo

ssave = 1;

tt = [tpl;tl]'; ii = [xpl(:,2);xl(:,2)]';

riepilogo(data,tt,ii,ssave);

if stat == 0
    return
else
    clear options
end

%% statistica

data(1).infSim = [tpl, xpl(:,2)];
data(2).infSim = [tl, xl(:,2);];

statistica(data)
