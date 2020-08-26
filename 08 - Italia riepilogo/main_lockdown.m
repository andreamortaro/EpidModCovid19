%{
Dinamica cambiata: il 09/03/20 parte il lockdown e nel sir entra in gioco
il fattore di controllo K.
Aggiustamento versione 2020-06-01: trovo un K per ogni giorno di lockdown e
ricavo un K(t) con t=giorno. Provo poi a fittare tale K(t) con una funzione
continua da inserire in un secondo momento nel modello.
%}

function [t,x,A] = main_lockdown(k0_c, pnt, kguess, deltatc, ffunz, oold, ffig, ssave)

%
%   [t,x,A] = main_lockdown(k0_c,kspan, window, kguess, pplot, oold, ffig, ssave))
%
%   Per ogni punto di kspan integro su una finestra, dove i parametri sono
%   descritti nella struttura window, e ricavo il valore del parametro di controllo k.
%   Trovo un k per ogni punto (=giorno) di kspan.
%
%   INPUTS:
%   kspan       : Intervallo temporale considerato.
%   window      : struttura che contiene i seguenti campi:
%                 - window.kl   : estremo sinistro finestra integrazione rispetto t_i
%                 - window.kr   : estremo destro finestra integrazione rispetto t_i
%                 - window.k0   : guess iniziale
%                 - window.pnt  : moltiplicatore per nodi integrazione
%   k0_c        : guess.
%
%   OUTPUTS:
%   iter       : numero di iterazioni fatte nel ciclo for
%                (corrisponde con length(kspan) - kr*h)
%   days       : il giorno days(i) corrisponde al k_c(i) parametro trovato
%   k_c        : parametri discreti trovati per ogni punto di kpsan
%
%

global t_u t_c Nass Ibar Rbar beta gamma date

% imposto latex come inteprete per i grafici
set(groot,...
    'defaulttextinterpreter','latex',...
    'defaultAxesTickLabelInterpreter','latex',...
    'defaultLegendInterpreter','latex');

%% PLOT FUNZIONALE

if ffunz == 1
    
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
    [~,z] = zoomPlot(K1,L1,[1e-12 1e-2],[0.4 0.255 0.4 0.4],[1 3]);
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


%% Lockdown: 1. stimo i K discreti

% definisco la finestra
window.h   = 1;              % daily time step
window.kl  = 3;
window.kr  = 4;
kspan      = t_u:1:t_c-window.kr*window.h;     % intervallo ricerca k discreto

% per comodità, procedura per i k_c discreti su kspan su un file a parte
[days, k_c] = stima_kdiscreti(kspan,window,k0_c,pnt);

T = table(days,k_c,'VariableNames',{'t_i' 'k_c(t_i)'}) %#ok<NOPRT>

%% PARTE 2 : Fitto i k_c discreti ottenuti e ricavo k(t)

problem2.options    = optimoptions('fmincon','Display','iter');
problem2.solver     = 'fmincon';
problem2.objective  = @minquad_kcontinuo;           % funzionale obiettivo minimizzare
problem2.x0         = kguess;                      % guess iniziale
%problem2.nonlcon = @(A)mycon(A);                  % vincolo non lineare su k (=beta>0)

% tecnicamente il vincolo non lineare dovrei imporlo, tuttavia
% l'ho imposto quando ho cercato i k discreti?

A = fmincon(problem2);

% update function
K = @(t) -A(1)*t.^2 + A(2)*t - A(3);


if ffig == 1
    fitting = figure();
    set(gca,'FontSize',12.5);
    nstep = 50;
    tt = linspace(t_u,t_c,nstep);
    plot(days,k_c,'*',tt,K(tt'));
    title("fitting $\kappa$");
    text(14,max(k_c)*.75,["a= " num2str(A(1)) "b= " num2str(A(2)) "c= " num2str(A(3)) ])
    xlabel("t (days)")
    ylabel("$\kappa$")

    if ssave == 1
        exportgraphics(fitting,'fittingk.pdf','ContentType','vector',...
                       'BackgroundColor','none')
    end
end
    
%% Simulazione modello oltre il lockdown

% update sistema
SI = @(t,x) [-(beta - x(1)*x(2)/K(t))*x(1)*x(2);
              (beta - x(1)*x(2)/K(t))*x(1)*x(2) - gamma*x(2)];

Jac = @(t,x) [ -beta*x(2) + 2*x(1)*(x(2)^2)/K(t), -beta*x(1) + 2*(x(1)^2)*x(2)/K(t);
                beta*x(2) - 2*x(1)*(x(2)^2)/K(t),  beta*x(1) - 2*(x(1)^2)*x(2)/K(t) - gamma];
options.Jacobian = Jac;

I0 = Ibar(t_u+1); R0 = Rbar(t_u+1); S0 = Nass-I0-R0;
x0 = [S0;I0]/Nass;          % dato iniziale in percentuale

nstep = 50;
[t, x]  = eulerorosenbrock(SI,linspace(t_u,t_c+deltatc,nstep),x0,options);
x(:,3) = ones(length(t),1) - x(:,1) - x(:,2);      % ricavo R per post-processing
x = Nass.*x;                                       % normalizzo

if ffig == 1
    
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
end

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