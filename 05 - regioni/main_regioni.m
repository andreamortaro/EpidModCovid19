close all
clear all
clc

% Identificazione parametri

global t_0 t_u t_c date regione Ibar Rbar Nass beta gamma days k_c

% Upload dati protezione civile
tmp = fullfile('..','00 - dpc_data','2020-05-22','dati-regioni');
[status,result] = fileattrib(tmp);
path_folder = result.Name;              % percorso alla cartella
[reg_label,date,Ibar_allreg,Rbar_allreg] = data_read_dpc_regioni(path_folder);
% Attenzione: poiche' gli indici in matrice partono da 1 e tm parte da 0
% date(i), Ibar(i), Rbar(i) e' in corrispondenza con t_i-1

N_reg = (1:21)';
T1 = table(N_reg,reg_label)

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

%% Pre-Lockdown [t_0,t_u]: stima beta e gamma prima del Lockdown
% (gestisco all'interno la figura)

K = stima_preLock();

beta = K(1); gamma = K(2);
R_0 = beta/gamma;
table(regione,beta,gamma,R_0)


%% Lockdown [t_u,t_c]: trovati beta e gamma calcolo il parametro di controllo k

% definisco la finestra con daily time step
h  = 1; kl = 3; kr = 4;
window.h = h;   window.kl = kl; window.kr = kr;

k0_c = 1;                  % guess iniziale (ottengo sempre gli stessi k_c)
kspan = t_u:1:t_c-kr*h;    % intervallo ricerca k discreto
options.pnt = 10;          % aumento numero nodi

% trovo i k_c discreti nell'intervallo [t_u, t_c]
[days, k_c] = stima_Lock(kspan,window,k0_c,options);

T2 = table(days,k_c,'VariableNames',{'t_i' 'k_c(t_i)'}) %#ok<NOPTS>


%% Fitting dei k_c discreti ottenuti

a = 1e-6;
b = 1e-4;
c = 1e-3;

A0 = [a,b,c];
A = fitting_k(A0,days,k_c);

%update K
K = @(t) -A(1)*t.^2 + A(2)*t - A(3);

%FIGURA
tspan1 = linspace(t_u,t_c,50);

fig_fittingk = figure();
hold on
plot(days,k_c,'*',tspan1,K(tspan1'));
box on
title({['fitting $\kappa$']
       ['a= ' num2str(A(1))]
       ['b= ' num2str(A(2))]
       ['c= ' num2str(A(3))]
       });
set(gca,'FontSize',12.5);

exportgraphics(fig_fittingk,'figure/' + regione + '_fittingk.pdf',...
               'ContentType','vector',...
               'BackgroundColor','none')

%% Simulazione modello con i parametri trovati (durante LOCKDOWN)

% update sistema
SI = @(t,x) [-(beta - x(1)*x(2)/K(t))*x(1)*x(2);
              (beta - x(1)*x(2)/K(t))*x(1)*x(2) - gamma*x(2)];

Jac = @(t,x) [ -beta*x(2) + 2*x(1)*(x(2)^2)/K(t), -beta*x(1) + 2*(x(1)^2)*x(2)/K(t);
                beta*x(2) - 2*x(1)*(x(2)^2)/K(t),  beta*x(1) - 2*(x(1)^2)*x(2)/K(t) - gamma];
options.Jacobian = Jac;

I0 = Ibar(t_u+1); R0 = Rbar(t_u+1); S0 = Nass-I0-R0;
x0 = [S0;I0]/Nass;          % dato iniziale in percentuale

nstep = 100;
tspan2 = linspace(t_u,t_c,nstep);            % durante il Lockdown

[t, x]  = eulerorosenbrock(SI,tspan2,x0,options);
x(:,3) = ones(length(t),1) - x(:,1) - x(:,2);      % ricavo R per post-processing
x = Nass.*x;      % normalizzo


fig_Lock = figure();
%plot(t,x(:,2),'r',t,x(:,3),'b',tm,Ibar(tm+1),'rx',tm,Rbar(tm+1),'bx');
plot(t,x(:,2),'SeriesIndex',1);
hold on
plot(t,x(:,3),'SeriesIndex',2);

tt = t_u:1:t_c;
plot(tt,Ibar(tt+1),'SeriesIndex',1,'LineStyle','none','Marker','*');
plot(tt,Rbar(tt+1),'SeriesIndex',2,'LineStyle','none','Marker','*');

ax = gca;
ax.XTick = [t_u,37,67,t_c];
ax.XTickLabel = date([t_u,37,67,t_c]+1);
ax.XTickLabelRotation = 45;

box on
legend('I','R','$I_{bar}$','$R_{bar}$','Location','NorthWest');
title([char(regione), ', simulazione SIR controllato']);
set(gca,'FontSize',12.5);

exportgraphics(fig_Lock,'figure/' + regione + '_Lock.pdf',...
               'ContentType','vector',...
               'BackgroundColor','none')
