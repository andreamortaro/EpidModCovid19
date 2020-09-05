function [K,A] = stima_kcontinuo(data,K0_cont,ffig,ssave)

%
%   [K,A] = stima_kcontinuo(data,K0_cont,ffig,ssave)
%
%   Fitto i parametri K_disc che si trovano in data.
%
%   INPUTS:
%   data        : struttura contenente i dati utili
%   K0_cont     : guess iniziale
%   ffig        : flag sulla stampa della figura
%   ssave       : flag print della figura
%
%   OUTPUTS:
%   K           : funzione risultante dall'ottimizzazione, fitta i K_disc
%   A           : parametri stimati dall'ottimizzazione
%

global days K_disc

% recupero i valori che servono
[~,~,t_u,t_c,date] = data.time;
[days,K_disc] = data.fittingK;

if isfield(data,'regione')
    regione = data(1).regione;
end

problem2.options    = optimoptions('fmincon','Display','iter');
problem2.solver     = 'fmincon';
problem2.objective  = @minquad_kcontinuo;          % funzionale obiettivo minimizzare
problem2.x0         = K0_cont;                     % guess iniziale

A = fmincon(problem2);

K = @(t) A(1)*exp(-((t-A(2))/A(3)).^2);  % guassiana

if ffig == 1
    
    % imposto latex come inteprete per i grafici
    set(groot,...
        'defaulttextinterpreter','latex',...
        'defaultAxesTickLabelInterpreter','latex',...
        'defaultLegendInterpreter','latex');

    fitting = figure();
    
    tt = linspace(t_u,t_c,50);
    plot(days,K_disc,'o',...
            'MarkerSize',4,...
            'MarkerEdgeColor',[.5 .7 .1],...
            'MarkerFaceColor',[.8 .9 0]);
    
    hold on
    p2 = plot(tt,K(tt'),'black','LineWidth',2.5);
    p2.Color(4) = 0.7;
    ylabel("$\kappa$")
    
    if exist('regione','var') == 1
        title(char(regione));
    else
        title('Italia')
    end

    ax = gca;
    ax.XTick = [t_u,37,67,98,t_c];
    ax.XTickLabel = date([t_u,37,67,98,t_c]+1);
    ax.XTickLabelRotation = 45;
    
    set(gca,'FontSize',12.5);

   
    if ssave == 1
        
        if exist('regione','var') == 1
            exportgraphics(fitting,'figure/' + regione + '/fittingk.pdf',...
            'ContentType','vector',...
            'BackgroundColor','none')
        else
            exportgraphics(fitting,'figure/fittingk.pdf',...
            'ContentType','vector',...
            'BackgroundColor','none')
        end
        
    end
end

end