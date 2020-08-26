% Identificazione parametri:

close all
clear
clc

global  x0 tm ym Nass t_0 t_u pnt %PRE LOCKDOWN

global Nass beta gamma t_u t_c Ibar Rbar days k_c   %LOCKDWON

ssave = 1;  % salvare le figure

% Upload dati protezione civile
tmp = fullfile('..','00 - dpc_data','2020-05-22','dati-andamento-nazionale');
%tmp = fullfile('..','00 - dpc_data','2020-06-30','dati-andamento-nazionale');

[status,result] = fileattrib(tmp);
path_folder = result.Name;              % percorso alla cartella
[date,Ibar,Rbar] = data_read_dpc(path_folder);
% Attenzione: poiche' gli indici in matrice partono da 1 e tm parte da 0
% date(i), Ibar(i), Rbar(i) e' in corrispondenza con t_i-1

% imposto latex come inteprete per i grafici
set(groot,...
    'defaulttextinterpreter','latex',...
    'defaultAxesTickLabelInterpreter','latex',...
    'defaultLegendInterpreter','latex');

% DATI:
t_0 = 0;                           % 2020-02-24
t_u = 14;                          % 2020-03-09 t finale senza controllo
t_c = length(date)-1;              % decremento perche' parto da 0
Nass = 60317000;                   % popolazione italiana istat 11.02.2020

%% Pre-Lockdown

tm  = t_0:1:t_u;                    % tm = [0,..,t_u]
ym  = [Ibar(tm+1),Rbar(tm+1)];

K0  = [0.5,0.1];                    % guess iniziale per [beta,gamma]
pnt = 5;                            % piu nodi per migliore risoluzione sistema minquad

I0 = Ibar(t_0+1); R0 = Rbar(t_0+1); S0 = Nass-I0-R0;
x0 = [S0;I0]/Nass;                  % dato iniziale in percentuale

% Minimizzazione

problem.options = optimoptions('fmincon','Display','iter');
problem.solver  = 'fmincon';
problem.objective = @minquad;
problem.x0 = K0;
problem.lb = [0,0];

K = fmincon(problem)    %#ok<NOPTS>
k1 = K(1); k2 = K(2);   % update parametri stimati 

SI = @(t,x) [-k1*x(1)*x(2); k1*x(1)*x(2) - k2*x(2)];   % update sys, nuovi parametri stimati
Jac = @(t,x) [-k1*x(2), -k1*x(2);
                            k1*x(2), k1*x(1) - k2];
options.Jacobian = Jac;

[t, x]  = rk4(SI,[t_0,t_u],x0,options);         % sistema lineare
x(:,3) = ones(length(t),1) - x(:,1) - x(:,2);   % ricavo R per post-processing
x = Nass.*x;                                    % ri-normalizzo da percentuale a Nass

% salvo i risultati
tpl = t;
xpl = x;

beta  = k1;
gamma = k2;

%% Lockdown: 1. stimo i K discreti

% definisco la finestra
h  = 1;                     % daily time step
kl = 3;
kr = 4;
window.h = h;   window.kl = kl; window.kr = kr;

k0_c = 1e-5;                % guess iniziale (ottengo sempre gli stessi k_c)
kspan = t_u:1:t_c-kr*h;     % intervallo ricerca k discreto
options.pnt = 100;          % aumento numero nodi integrazione

% per comodità, procedura per i k_c discreti su kspan su un file a parte
[days, k_c] = stima_kdiscreti(kspan,window,k0_c,options);

%% Lockdown: 2. Fitto i k_c discreti ottenuti e ricavo k(t)

a = 1e-6; b = 1e-4; c = 1e-3;

problem2.options    = optimoptions('fmincon','Display','iter');
problem2.solver     = 'fmincon';
problem2.objective  = @minquad_kcontinuo;           % funzionale obiettivo minimizzare
problem2.x0         = [a,b,c];                      % guess iniziale
%problem2.nonlcon = @(A)mycon(A);                  % vincolo non lineare su k (=beta>0)

A = fmincon(problem2);

% update function
K = @(t) -A(1)*t.^2 + A(2)*t - A(3);

%% Simulazione modello oltre il lockdown

% update sistema
SI = @(t,x) [-(beta - x(1)*x(2)/K(t))*x(1)*x(2);
              (beta - x(1)*x(2)/K(t))*x(1)*x(2) - gamma*x(2)];

Jac = @(t,x) [ -beta*x(2) + 2*x(1)*(x(2)^2)/K(t), -beta*x(1) + 2*(x(1)^2)*x(2)/K(t);
                beta*x(2) - 2*x(1)*(x(2)^2)/K(t),  beta*x(1) - 2*(x(1)^2)*x(2)/K(t) - gamma];
options.Jacobian = Jac;

I0 = Ibar(t_u+1); R0 = Rbar(t_u+1); S0 = Nass-I0-R0;
x0 = [S0;I0]/Nass;          % dato iniziale in percentuale

tspan = linspace(t_u,t_c+10,50);
[t, x]  = eulerorosenbrock(SI,tspan,x0,options);
x(:,3) = ones(length(t),1) - x(:,1) - x(:,2);      % ricavo R per post-processing
x = Nass.*x;      % normalizzo

% salvo i risultati
tl = t;
xl = x;

%% Figure: riepilogo

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
%axis([0 3e5])      % come albi

box on
legend('I','$I_{bar}$','Location','Best');
ylabel('casi confermati');
set(gca,'FontSize',12.5)


if ssave == 1
    exportgraphics(fig,['italia_riepilogo ' num2str(date(t_c+1)) '.pdf'],'ContentType','vector',...
                   'BackgroundColor','none')
end

% VINCOLO NON LINEARE SU K Quando voglio fittarei k discreti con una
% function k(t) continua

function [c,ceq] = mycon(A)

global x0 beta gamma t_u t_c Ibar Rbar Nass

    K = @(t) A(1)*exp(A(2)*t).*(1-exp(A(3)*t));
 
    SI = @(t,x) [-(beta - x(1)*x(2)/K(t))*x(1)*x(2);
                  (beta - x(1)*x(2)/K(t))*x(1)*x(2) - gamma*x(2)];
          
    Jac = @(t,x) [ -beta*x(2) + 2*x(1)*(x(2)^2)/K(t), -beta*x(1) + 2*(x(1)^2)*x(2)/K(t);
                    beta*x(2) - 2*x(1)*(x(2)^2)/K(t),  beta*x(1) - 2*(x(1)^2)*x(2)/K(t) - gamma];
    options.Jacobian = Jac;
    
    I0 = Ibar(t_u+1); R0 = Rbar(t_u+1); S0 = Nass-I0-R0;
    x0 = [S0;I0]/Nass;          % dato iniziale in percentuale
    
    tspan = linspace(t_u,t_c,74);            % cosi non ho problemi con matrici singolari
    [t, xm]  = eulerorosenbrock(SI,tspan,x0,options);

    %xm = Nass.*xm;                                        % ???
    
    % controllo con valori in percentuale, altrimenti impossibile verificare la condizione
    
	% Nonlinear inequality constraints (c(K)<=0)
    
    if t ~= tspan'
        warning('attenzione ai vincoli')
    end
    
    c = xm(:,1).*xm(:,2) - 2*beta*K(t);   % deve essere <=0 (ATTENZIONE all'=)
    ceq = [];
end
