function [t, x, beta, gamma] = prelock(data, K0, options)

%
%   [t, x, beta, gamma] = prelock(data, K0, options)
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
%   beta        : parametro del SIR ottenuto dall'ottimizzazione
%   gamma       : parametro del SIR ottenuto dall'ottimizzazione
%

global t_0 t_u Nass Ibar Rbar date x0 tm ym pnt

% recupero i valori che servono
[Nass,Ibar,Rbar] = data.value;
[t_0,t_u,~,date] = data.time;

if nargin == 2
    ssens = 1;
    ffig = 1;
    ssave = 1;
else
    if isfield(options,'ssens')
        ssens = options.ssens;
    end
    if isfield(options,'ffig')
        ffig = options.ffig;
    end
    if isfield(options,'ssave')
        ssave = options.ssave;
    end 
end

% punti dove minimizzare
tm  = t_0:1:t_u;                    % tm = [0,..,t_u]
ym  = [Ibar(tm+1),Rbar(tm+1)];

I0 = Ibar(t_0+1); R0 = Rbar(t_0+1); S0 = Nass-I0-R0;
x0 = [S0;I0]/Nass;                  % dato iniziale in percentuale

pnt = 1;    % default
if(nargin==3)
    pnt = options.pnt;
end

% Minimizzazione

problem.options     = optimoptions('fmincon','Display','iter');
problem.solver      = 'fmincon';
problem.objective   = @minquad_prelock;
problem.x0          = K0;
problem.lb          = [0,0];

K = fmincon(problem) %#ok<NOPRT>
beta = K(1); gamma = K(2);   % update parametri stimati 

SI = @(t,x) [-beta*x(1)*x(2); beta*x(1)*x(2) - gamma*x(2)];   % update sys, nuovi parametri stimati

options.Jacobian = @(t,x) [-beta*x(2), -beta*x(2);
                            beta*x(2), beta*x(1) - gamma];

[t, x]  = rk4(SI,[t_0,t_u],x0,options);     % sistema lineare
x(:,3) = ones(length(t),1) - x(:,1) - x(:,2);   % ricavo R per post-processing
x = Nass.*x;                                    % ri-normalizzo da percentuale a Nass


%% FIGURA

if ffig == 1
    
    % imposto latex come inteprete per i grafici
    set(groot,...
        'defaulttextinterpreter','latex',...
        'defaultAxesTickLabelInterpreter','latex',...
        'defaultLegendInterpreter','latex');
    
    fig = figure();
    p1 = plot(tm,Ibar(tm+1),'o',...
        'MarkerSize',4,...
        'MarkerEdgeColor','red',...
        'MarkerFaceColor',[1 .6 .6]);        

    hold on
    p2 = plot(tm,Rbar(tm+1),'o',...
        'MarkerSize',4,...
        'MarkerEdgeColor',[.3 .4 .6],...
        'MarkerFaceColor',[.3 .6 .8]);
    
    p3 = plot(t,x(:,2),'SeriesIndex',2,'Linewidth',2.5);
    %p3 = plot(t,x(:,2),'color','black','Linewidth',2.5);

    hold on
    p4 = plot(t,x(:,3),'SeriesIndex',1,'Linewidth',2.5);
    
    p3.Color(4) = 0.6;
    p4.Color(4) = 0.6;
    
    ax = gca;
    ax.XTick = 0:7:14;
    ax.XTickLabel = date((0:7:14)+1);
    ax.XTickLabelRotation = 45;
    box on
    legend([p1,p2,p3,p4],'$I_{bar}$','$R_{bar}$','I','R','Location','NorthWest');
    ylabel('casi confermati');
    set(gca,'FontSize',12.5)
    
    if ssave == 1
        exportgraphics(fig,'figure/italia-preLock.pdf','ContentType','vector',...
                       'BackgroundColor','none')
    end
end

%% Analisi sensitivita'

if ssens == 1
    x1 = x(:,1)/Nass;	% S(t)
    x2 = x(:,2)/Nass;  	% I(t)
    x3 = x(:,3)/Nass;  	% R(t)

    D0 = [0 0 0 0 0 0]';   % condizione iniziale del sistema D

    % Matrice sensitivita'
    % m non e' il tempo, ma l'indice che gira sul vettore x

    sensys = @(m,D) [-beta*x2(m)*D(1)-beta*x1(m)*D(3)-x1(m)*x2(m);...
                     -beta*x2(m)*D(2)-beta*x1(m)*D(4);...
                      beta*x2(m)*D(1) + (beta*x1(m)-gamma)*D(3) + x1(m)*x2(m);...
                      beta*x2(m)*D(2) + (beta*x1(m)-gamma)*D(4) - x2(m);...
                      gamma*D(3);...
                      gamma*D(4) + x2(m)];

    % Risolvo con Eulero Esplicito, il sistemare e' lineare in D
    D = odesol(sensys,t,D0);

    fig2 = figure();
    plot(t,D(:,1),t,D(:,2),t,D(:,3),t,D(:,4),t,D(:,5),t,D(:,6))
    %title('Sensitivity')
    H=legend('$\partial_{\beta}S$','$\partial_{\gamma}S$',...
             '$\partial_{\beta}I$','$\partial_{\gamma}I$',...
             '$\partial_{\beta}R$','$\partial_{\gamma}R$');
    set(H,'Location','Best',...
          'NumColumns',2,...
          'Orientation','horizontal')
    axis tight
    grid on
    xlabel('t')
    set(gca,'FontSize',12.5)
    
    if ssave == 1
        exportgraphics(fig2,'figure/sensitivita.pdf','ContentType','vector',...
                       'BackgroundColor','none')
    end
               
end

end