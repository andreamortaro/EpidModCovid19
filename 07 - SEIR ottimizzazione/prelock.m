function [t, x, beta, gamma,mu] = prelock(data, K0, options)

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

% recupero i valori che servono
[Nass,Ibar,Rbar,~,Ebar] = data.value;
[t_0,t_1,t_u,~,date] = data.time;

if nargin == 2
    ffig = 1;
    ssave = 1;
    prelockopt =1;
else
    if isfield(options,'ffig')
        ffig = options.ffig;
    end
    if isfield(options,'ssave')
        ssave = options.ssave;
    end
end

if isfield(data,'regione')
    regione = data(1).regione;
    ssens = 0;      % nelle regioni non ho previsto analisi sensitivita
else    % in questo caso son sicuro di essere in italia
    if isfield(options,'ssens')
        ssens = options.ssens;
    else
        ssens = 1;
    end
end

pnt = 1;    % default
if(nargin==3)
    pnt = options.pnt;
end

%% ottimizzazione per beta e gamma

data(1).prelock = t_0;  % tstart
data(2).prelock = t_u;  % tfinal

[~,~,beta,gamma,mu] = stima_beta_gamma(data,K0,pnt);
        
%% simulazione

% update sys, nuovi parametri trovati
SEI = @(t,x) [-beta*x(1)*x(3);...
               beta*x(1)*x(3) - mu*x(2);...
               mu*x(2)- gamma*x(3)];

opt.InitialStep = 0.01;

E0 = Ebar(t_0+1); I0 = Ibar(t_0+1); R0 = Rbar(t_0+1); S0 = Nass-E0-I0-R0;
x0 = [S0;E0;I0]/Nass;                           % dato iniziale in percentuale

[t, x]  = rk4(SEI,[t_0,t_u],x0,opt);       % simulazione modello

x(:,4) = ones(length(t),1) - x(:,1) - x(:,2) - x(:,3);   % ricavo R per post-processing
x = Nass.*x;                                    % ri-normalizzo da percentuale a Nass

tm = t_0:1:t_u;

if ffig == 1
    
    % imposto latex come inteprete per i grafici
    set(groot,...
        'defaulttextinterpreter','latex',...
        'defaultAxesTickLabelInterpreter','latex',...
        'defaultLegendInterpreter','latex');
    
    fig = figure();
    p1 = plot(tm,Ebar(tm+1),'o',...
        'MarkerSize',4,...
        'MarkerEdgeColor','green',...
        'MarkerFaceColor',[1 .6 .6]);        

    hold on
    p2 = plot(tm,Ibar(tm+1),'o',...
        'MarkerSize',4,...
        'MarkerEdgeColor','red',...
        'MarkerFaceColor',[1 .6 .6]); 
    
    p3 = plot(tm,Rbar(tm+1),'o',...
        'MarkerSize',4,...
        'MarkerEdgeColor',[.3 .4 .6],...
        'MarkerFaceColor',[.3 .6 .8]);
    
    p4 = plot(t,x(:,2),'SeriesIndex',3,'Linewidth',2.5); p3.Color(4) = 0.6;
    p5 = plot(t,x(:,3),'color','black','Linewidth',2.5); p4.Color(4) = 0.6;
    p6 = plot(t,x(:,4),'SeriesIndex',1,'Linewidth',2.5); p4.Color(4) = 0.6;

    ax = gca;
    ax.XTick = 0:7:14;
    ax.XTickLabel = date((0:7:14)+1);
    ax.XTickLabelRotation = 45;
    box on
    legend([p1,p2,p3,p4,p5,p6],'$E_{bar}$','$I_{bar}$','$R_{bar}$',...
                                'E','I','R','Location','NorthWest');
    ylabel('casi attuali confermati');
    
    if exist('regione','var') == 1
        title(char(regione));
        limsy=get(gca,'YLim');
        set(gca,'Ylim',[0 limsy(2)]);
    else
        title('Italia')
    end
    
    if ssave == 1
        if exist('regione','var') == 1
            exportgraphics(fig,'figure/' + regione + '/prelock.pdf','ContentType','vector',...
            'BackgroundColor','none')
        else
            exportgraphics(fig,'figure/Italia/prelock.pdf','ContentType','vector',...
            'BackgroundColor','none')
        end
    end
end

% %% Analisi sensitivita'
% 
% if ssens == 1
%     x1 = x(:,1)/Nass;	% S(t)
%     x2 = x(:,2)/Nass;  	% I(t)
%     x3 = x(:,3)/Nass;  	%#ok<NASGU> % R(t)
% 
%     D0 = [0 0 0 0 0 0]';   % condizione iniziale del sistema D
% 
%     % Matrice sensitivita'
%     % m non e' il tempo, ma l'indice che gira sul vettore x
% 
%     sensys = @(m,D) [-beta*x2(m)*D(1)-beta*x1(m)*D(3)-x1(m)*x2(m);...
%                      -beta*x2(m)*D(2)-beta*x1(m)*D(4);...
%                       beta*x2(m)*D(1) + (beta*x1(m)-gamma)*D(3) + x1(m)*x2(m);...
%                       beta*x2(m)*D(2) + (beta*x1(m)-gamma)*D(4) - x2(m);...
%                       gamma*D(3);...
%                       gamma*D(4) + x2(m)];
% 
%     % Risolvo con Eulero Esplicito, il sistemare e' lineare in D
%     D = odesol(sensys,t,D0);
% 
%     fig2 = figure();
%     p1 = plot(t,D(:,1),'Linewidth',1.5);
%     hold on
%     p2 = plot(t,D(:,2),'Linewidth',1.5);
%     p3 = plot(t,D(:,3),'Linewidth',1.5);
%     p4 = plot(t,D(:,4),'Linewidth',1.5);
%     p5 = plot(t,D(:,5),'Linewidth',1.5);
%     p6 = plot(t,D(:,6),'Linewidth',1.5);
%     %title('Sensitivity')
%     H=legend([p1,p2,p3,p4,p5,p6],'$\partial_{\beta}S$','$\partial_{\gamma}S$',...
%              '$\partial_{\beta}I$','$\partial_{\gamma}I$',...
%              '$\partial_{\beta}R$','$\partial_{\gamma}R$');
%     set(H,'Location','Best',...
%           'NumColumns',2,...
%           'Orientation','horizontal')
%     axis tight
%     grid on
%     xlabel('t')
%     
%     set(gca,'FontSize',12.5)
%     if ssave == 1
%         exportgraphics(fig2,'figure/Italia/sensitivita.pdf','ContentType','vector',...
%                        'BackgroundColor','none')
%     end
%                
% end
% 
% end
