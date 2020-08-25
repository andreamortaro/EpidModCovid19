%{
Dinamica cambiata: il 09/03/20 parte il lockdown e nel sir entra in gioco
il fattore di controllo K.
Aggiustamento versione 2020-06-01: trovo un K per ogni giorno di lockdown e
ricavo un K(t) con t=giorno. Provo poi a fittare tale K(t) con una funzione
continua da inserire in un secondo momento nel modello.
%}

%% Identificazione parametri:

close all
clear all
clc

global Nass beta gamma t_u t_c Ibar Rbar days k_c

pplot = 1;  % stampare funzionale
oold  = 1;  % stampo il funzionale nella vecchia maniera
ssave = 1;  % salvare le figure


% Upload dati protezione civile
tmp = fullfile('..','00 - dpc_data','2020-05-22','dati-andamento-nazionale');
[status,result] = fileattrib(tmp);
path_folder = result.Name;              % percorso alla cartella
[date,Ibar,Rbar] = data_read_dpc(path_folder);
% Attenzione: poiche' gli indici in matrice partono da 1 e tm parte da 0
% date(i), Ibar(i), Rbar(i) e' in corrispondenza con t_i-1

% DATI:
t_u = 14;                          % 2020-03-09: inizio lockdown
t_c = length(date)-1;              % decremento perche' parto da 0

% parametri trovati prima con main.m
beta    = 0.30883;
gamma   = 0.04953;
Nass    = 60317000;             % popolazione italiana istat 11.02.2020
%Nass = 6e9;                    % imbrogliando così arrivo a k~1e-5

% imposto latex come inteprete per i grafici
set(groot,...
    'defaulttextinterpreter','latex',...
    'defaultAxesTickLabelInterpreter','latex',...
    'defaultLegendInterpreter','latex');

%% PLOT FUNZIONALE

if pplot == 1
    
    if oold == 1
        % Primo intevallo
        kspan = linspace(1e-12,1e-2,50);
        pnt = 10;
        [K1,L1] = plot_minquad_kdiscreti(kspan,pnt);

        fig1 = figure();
        set(gca,'FontSize',12.5);
        plot(K1,L1)
        xlabel("$\kappa$")
        ylabel("L($\kappa$)")
        title(['Funzionale in [' num2str(kspan(1)) ',' num2str(kspan(end)) ']'])
        grid on


        if ssave == 1
            exportgraphics(fig1,'funzionale1.pdf','ContentType','vector',...
                           'BackgroundColor','none')
        end
        
        % secondo plot  
        kspan = linspace(1e-2,10,150);
        pnt = 15;
        [K2,L2] = plot_minquad_kdiscreti(kspan,pnt);

        fig2 = figure();
        set(gca,'FontSize',12.5);
        plot(K2,L2)
        xlabel("$\kappa$")
        ylabel("L($\kappa$)")
        title(['Funzionale in [' num2str(kspan(1)) ',' num2str(kspan(end)) ']'])
        grid on
        hold on

        if ssave == 1
            exportgraphics(fig2,'funzionale2.pdf','ContentType','vector',...
                           'BackgroundColor','none')
        end
    end
    
    % plot con zoom
    kspan = linspace(1e-12,0.05,150);
    pnt = 10;
    [K1,L1] = plot_minquad_kdiscreti(kspan,pnt);

    fig = figure();
    xlabel("$\kappa$")
    ylabel("L($\kappa$)")
    set(gca,'FontSize',12.5);
    grid on
    box on
    hold on 
    
    % riplotto sopra i risultati ma su un intervallo diverso
    kspan = linspace(1e-12,10,100);
    pnt = 15;
    [K2,L2] = plot_minquad_kdiscreti(kspan,pnt);
    plot(K2,L2)
    
    % zoomPlot to highlight a portion of the major plot
    [p,z] = zoomPlot(K1,L1,[1e-12 1e-2],[0.4 0.255 0.4 0.4],[1 3]);
    z.XLim = [0 0.01];
    z.YLim = [0 1e14];
    z.XGrid = 'on'; z.YGrid = 'on';
    z.XTick = 0:0.0025:0.010;
	z.XAxis.TickLabelFormat = '%,.3f';
    set(z,'FontSize',10);
    
    hold on
    legend hide
    
    if ssave == 1
        exportgraphics(fig,'funzionalezoom.pdf','ContentType','vector',...
                       'BackgroundColor','none')
    end      
end


%% Stimo i K discreti

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

T = table(days,k_c,'VariableNames',{'t_i' 'k_c(t_i)'}) %#ok<NOPTS>


%% PARTE 2 : Fitto i k_c discreti ottenuti e ricavo k(t)

a = 1e-6;
b = 1e-4;
c = 1e-3;

% a =   3.387e-06;  %(3.319e-06, 3.456e-06)
% b =   0.0003582;  %(0.0003515, 0.000365)
% c =    0.002778;  %(0.00263, 0.002926)


problem2.options    = optimoptions('fmincon','Display','iter');
problem2.solver     = 'fmincon';
problem2.objective  = @minquad_kcontinuo;           % funzionale obiettivo minimizzare
problem2.x0         = [a,b,c];                      % guess iniziale
%problem2.nonlcon = @(A)mycon(A);                  % vincolo non lineare su k (=beta>0)

% tecnicamente il vincolo non lineare dovrei imporlo, tuttavia
% l'ho imposto quando ho cercato i k discreti?

A = fmincon(problem2);

% update function
K = @(t) -A(1)*t.^2 + A(2)*t - A(3);

tspan = linspace(t_u,t_c,50);

figure()
set(gca,'FontSize',12.5);
plot(days,k_c,'*',tspan,K(tspan'));
title("fitting $\kappa$");
text(14,max(k_c)*.75,["a= " num2str(A(1)) "b= " num2str(A(2)) "c= " num2str(A(3)) ])
xlabel("t (days)")
ylabel("$\kappa$")

if ssave == 1
    exportgraphics(fig,'fittingk.pdf','ContentType','vector',...
                   'BackgroundColor','none')
end

% FIGURA

% update sistema
SI = @(t,x) [-(beta - x(1)*x(2)/K(t))*x(1)*x(2);
              (beta - x(1)*x(2)/K(t))*x(1)*x(2) - gamma*x(2)];

Jac = @(t,x) [ -beta*x(2) + 2*x(1)*(x(2)^2)/K(t), -beta*x(1) + 2*(x(1)^2)*x(2)/K(t);
                beta*x(2) - 2*x(1)*(x(2)^2)/K(t),  beta*x(1) - 2*(x(1)^2)*x(2)/K(t) - gamma];
options.Jacobian = Jac;

I0 = Ibar(t_u+1); R0 = Rbar(t_u+1); S0 = Nass-I0-R0;
x0 = [S0;I0]/Nass;          % dato iniziale in percentuale

[t, x]  = eulerorosenbrock(SI,tspan,x0,options);
x(:,3) = ones(length(t),1) - x(:,1) - x(:,2);      % ricavo R per post-processing
x = Nass.*x;      % normalizzo

tt = t_u:1:t_c;

fig = figure();
plot(t,x(:,2),'SeriesIndex',1,'LineWidth',2);
hold on
plot(t,x(:,3),'SeriesIndex',2,'LineWidth',2);
plot(tt,Ibar(tt+1),'o',...
    'MarkerSize',3,...
    'MarkerEdgeColor','blue',...
    'MarkerFaceColor',[1 .6 .6]);
plot(tt,Rbar(tt+1),'o',...
    'MarkerSize',3,...
    'MarkerEdgeColor','red',...
    'MarkerFaceColor',[1 .6 .6]);

ax = gca;
ax.XTick = [t_u,37,67,t_c];
ax.XTickLabel = date([t_u,37,67,t_c]+1);
ax.XTickLabelRotation = 45;

box on
legend('I','R','$I_{bar}$','$R_{bar}$','Location','NorthWest');
title("simulazione SIR controllato");

if ssave == 1
    exportgraphics(fig,'italia-lockdown.pdf','ContentType','vector',...
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
