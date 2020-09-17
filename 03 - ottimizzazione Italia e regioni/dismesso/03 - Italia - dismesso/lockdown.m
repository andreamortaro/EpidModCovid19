function [t,x,days,K_disc,Kfun,A] = lockdown(data, K0_disc, K0_cont, options)

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
%   K_disc      : k discreti trovati
%   Kfun        : funzione che fitta i K_disc
%   A           : parametri stimati per fittare i K_disc
%

% recupero i valori che servono
[Nass,Ibar,Rbar] = data.value;
[~,~,t_u,t_c,date] = data.time;
[beta,gamma] = data.parameters;

if nargin == 2
    ffunz = 1;
    oold = 1; %#ok<NASGU>
    ffig = 1;
    ssave = 1;
else
    if isfield(options,'ffunz')
        ffunz = options.ffunz;
    end
    if isfield(options,'oold')
        oold = options.oold; %#ok<NASGU>
    end    
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

%% 2. Fitto i K_disc discreti ottenuti e ricavo k(t) continua

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

nstep = 50;
deltatc = 0;
if isfield(options,'deltatc')
    deltatc = options.deltatc;
end

[t, x]  = eulerorosenbrock(SI,linspace(t_u,t_c+deltatc,nstep),x0,options);
x(:,3) = ones(length(t),1) - x(:,1) - x(:,2);      % ricavo R per post-processing
x = Nass.*x;                                       % normalizzo

if ffig == 1
    
    % imposto latex come inteprete per i grafici
    set(groot,...
        'defaulttextinterpreter','latex',...
        'defaultAxesTickLabelInterpreter','latex',...
        'defaultLegendInterpreter','latex');

    tt = t_u:1:t_c;
    
    fig = figure();
    ax_fig = axes;
    
    p1 = plot(tt,Ibar(tt+1),'o',...
        'MarkerSize',4,...
        'MarkerEdgeColor','red',...
        'MarkerFaceColor',[1 .6 .6]);

    hold on
    p2 = plot(tt,Rbar(tt+1),'o',...
        'MarkerSize',4,...
        'MarkerEdgeColor',[.3 .4 .6],...
        'MarkerFaceColor',[.3 .6 .8]);

    p3 = plot(t,x(:,2),'SeriesIndex',2,'LineWidth',2.5);

    p4 = plot(t,x(:,3),'SeriesIndex',1,'LineWidth',2.5);    
    p3.Color(4) = 0.6;
    p4.Color(4) = 0.6;
    
    ax = gca;
    ax.XTick = [t_u,37,67,t_c];
    ax.XTickLabel = date([t_u,37,67,t_c]+1);
    ax.XTickLabelRotation = 45;

    box on
    legend([p1,p2,p3,p4],'$I_{bar}$','$R_{bar}$','I','R','Location','NorthWest');
    ylabel('casi confermati');
    title('Italia');
        
    fig_xlim = get(gca,'XLim');
    
	set(gca,'FontSize',12.5);

    ax_fig.XLim = fig_xlim;
    ax_fig.YLim = [0 3.5e5];
    
    p = get(gca, 'Position');
    h = axes('Parent',gcf,'Position', [p(1)+.46 p(2)+.46 p(3)-.5 p(4)-.5],'box','on');

    hold(h,'on')
    
    %%% INSERISCO LA CURVA K
    
    % imposto latex come inteprete per i grafici
    set(groot,...
        'defaulttextinterpreter','latex',...
        'defaultAxesTickLabelInterpreter','latex',...
        'defaultLegendInterpreter','latex');

    set(gca,'FontSize',8);
    
    tt = linspace(t_u,t_c,50);
    plot(days,K_disc,'o',...
            'MarkerSize',3,...
            'MarkerEdgeColor',[.5 .7 .1],...
            'MarkerFaceColor',[.8 .9 0]);
    
    hold on
    p2 = plot(tt,Kfun(tt'),'black','LineWidth',2);
    p2.Color(4) = 0.7;
    ylabel("$\kappa$")

    ax = gca;
    ax.XTick = [t_u,37,67,t_c];
    ax.XTickLabel = date([t_u,37,67,t_c]+1);
    ax.XTickLabelRotation = 45;
        
    if ssave == 1
        exportgraphics(fig,['figure/italia_lockdown' num2str(date(t_c+1)) '.pdf'],...
                       'BackgroundColor','none');
    end
    
end

end
