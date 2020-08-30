function [t,x,A,days,K_disc] = regioni_lockdown(data, K0_disc, K0_cont, options)

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
%   A           : parametri stimati per fittare i K_disc
%   days        : giorni (non indici) in cui ho trovato k discreto
%   K_disc      : il valore di k discreto
%

global t_u t_c Nass Ibar Rbar beta gamma date K_disc days regione

% recupero i valori che servono
[Nass,Ibar,Rbar] = data.value;
[~,t_u,t_c,date] = data.time;
regione = data(1).regione;
[beta,gamma] = data.parameters;

if nargin == 2
    ffig = 1;
    ssave = 1;
else
    if isfield(options,'ffig')
        ffig = options.ffig;
    end
    if isfield(options,'ssave')
        ssave = options.ssave;
    end
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
[days, K_disc] = stima_kdiscreti(kspan,window,K0_disc,pnt);

%% 2. Fitto i K_disc discreti ottenuti e ricavo k(t)

problem2.options    = optimoptions('fmincon','Display','iter');
problem2.solver     = 'fmincon';
problem2.objective  = @minquad_kcontinuo;           % funzionale obiettivo minimizzare
problem2.x0         = K0_cont;                      % guess iniziale
problem2.lb         = [0,0,0];
%problem2.nonlcon    = @(A)mycon(A);                 % vincolo non lineare su k (=beta>0)
problem2.OptimalityTolerance = 1e-12;
problem2.StepTolerance = 1e-12;
problem2.FunctionTolerance = 1e-12;
% non cambia se non impongo mycon perche' gia' imposto per k disc

A = fmincon(problem2);

% update function
switch regione
    case {'Veneto','Emilia-Romagna','Piemonte'}
        Kfun = @(t) A(1)*exp(-((t-A(2))/A(3)).^2);
    otherwise
        Kfun = @(t) -A(1)*t.^2 + A(2)*t - A(3);
end

%FIGURA
if ffig == 1
    
    % imposto latex come inteprete per i grafici
    set(groot,...
        'defaulttextinterpreter','latex',...
        'defaultAxesTickLabelInterpreter','latex',...
        'defaultLegendInterpreter','latex');

    fitting = figure();
    set(gca,'FontSize',12.5);
    
    tt = linspace(t_u,t_c,50);
    plot(days,K_disc,'o',...
            'MarkerSize',4,...
            'MarkerEdgeColor',[.5 .7 .1],...
            'MarkerFaceColor',[.8 .9 0]);
    
    hold on
    p2 = plot(tt,Kfun(tt'),'black','LineWidth',2.5);
    p2.Color(4) = 0.7;
    %text(14,max(K_disc)*.75,["a= " num2str(A(1)) "b= " num2str(A(2)) "c= " num2str(A(3)) ])
    ylabel("$\kappa$")

    ax = gca;
    ax.XTick = [t_u,37,67,t_c];
    ax.XTickLabel = date([t_u,37,67,t_c]+1);
    ax.XTickLabelRotation = 45;
    
    title(char(regione),'FontSize',13);

    if ssave == 1
        exportgraphics(fitting,'figure/' + regione + '_fittingk.pdf',...
                       'ContentType','vector',...
                       'BackgroundColor','none')
    end
end
    
%% Simulazione modello oltre il lockdown

% update sistema
SI = @(t,x) [-(beta - x(1)*x(2)/Kfun(t))*x(1)*x(2);
              (beta - x(1)*x(2)/Kfun(t))*x(1)*x(2) - gamma*x(2)];

Jac = @(t,x) [ -beta*x(2) + 2*x(1)*(x(2)^2)/Kfun(t), -beta*x(1) + 2*(x(1)^2)*x(2)/Kfun(t);
                beta*x(2) - 2*x(1)*(x(2)^2)/Kfun(t),  beta*x(1) - 2*(x(1)^2)*x(2)/Kfun(t) - gamma];
options.Jacobian = Jac;

I0 = Ibar(t_u+1); R0 = Rbar(t_u+1); S0 = Nass-I0-R0;
x0 = [S0;I0]/Nass;          % dato iniziale in percentuale

deltatc = 0;
nstep = 500;
if isfield(options,'deltatc')
    deltatc = options.deltatc;
    nstep = options.nstep;
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
    %p3 = plot(t,x(:,2),'color','black','LineWidth',2.5);

    p4 = plot(t,x(:,3),'SeriesIndex',1,'LineWidth',2.5);    
    p3.Color(4) = 0.6;
    p4.Color(4) = 0.6;
    
    ax = gca;
    ax.XTick = [t_u,37,67,t_c];
    ax.XTickLabel = date([t_u,37,67,t_c]+1);
    ax.XTickLabelRotation = 45;

    box on
    legend([p1,p2,p3,p4],'$I_{bar}$','$R_{bar}$','I','R','Location','NorthWest');

    title(char(regione),'FontSize',13);

    if ssave == 1
        exportgraphics(fig,['figure/' char(regione) '_lockdown.pdf'],...
                       'ContentType','vector',...
                       'BackgroundColor','none')
    end
end

end

% VINCOLO NON LINEARE SU K Quando voglio fittarei k discreti con una
% function k(t) continua

function [c,ceq] = mycon(A)

global x0 beta gamma t_u t_c Ibar Rbar Nass regione

    switch regione
        case {'Veneto','Emilia-Romagna','Piemonte'}
            Kfun = @(t) A(1)*exp(-((t-A(2))/A(3)).^2);
        otherwise
            Kfun = @(t) -A(1)*t.^2 + A(2)*t - A(3);
    end
    
    SI = @(t,x) [-(beta - x(1)*x(2)/Kfun(t))*x(1)*x(2);
                  (beta - x(1)*x(2)/Kfun(t))*x(1)*x(2) - gamma*x(2)];

    Jac = @(t,x) [ -beta*x(2) + 2*x(1)*(x(2)^2)/Kfun(t), ...
                   -beta*x(1) + 2*(x(1)^2)*x(2)/Kfun(t);
                    beta*x(2) - 2*x(1)*(x(2)^2)/Kfun(t),...
                    beta*x(1) - 2*(x(1)^2)*x(2)/Kfun(t) - gamma];
    options.Jacobian = Jac;
    
    I0 = Ibar(t_u+1); R0 = Rbar(t_u+1); S0 = Nass-I0-R0;
    x0 = [S0;I0]/Nass;          % dato iniziale in percentuale
    
    tspan = linspace(t_u,t_c,30);            % cosi non ho problemi con matrici singolari
    [t, xm]  = eulerorosenbrock(SI,tspan,x0,options);

    %xm = Nass.*xm;                                        % ???
    
    % controllo con valori in percentuale, altrimenti impossibile verificare la condizione
    
	% Nonlinear inequality constraints (c(K)<=0)
    
    if t ~= tspan'
        warning('attenzione ai vincoli')
    end
    
    c = xm(:,1).*xm(:,2) - 2*beta*Kfun(t);   % deve essere <=0 (ATTENZIONE all'=)
    ceq = [];
end