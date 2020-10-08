function riepilogo(data,tt,ii,options)

%
%   la funzione riepilogo stampa l'intera simulazione prelock-lockdown del
%   nostro modello k-SIR.
%

[t_0,~,t_u,t_c,date] = data.time;
[~,Ibar,~] = data.value;
[days,K_disc,Kfun] = data.Kvalue;

ssave = 1;
if isfield(options,'ssave')
    ssave = options.ssave;
end

if isfield(data,'regione')
    regione = data(1).regione;
end

fig = figure();
ax_fig = axes;

% imposto latex come inteprete per i grafici
set(groot,...
    'defaulttextinterpreter','latex',...
    'defaultAxesTickLabelInterpreter','latex',...
    'defaultLegendInterpreter','latex');

p1 = plot(t_0:t_c,Ibar,'o',...
    'MarkerSize',4,...
    'MarkerEdgeColor','red',...
    'MarkerFaceColor',[1 .6 .6]);

hold on
p2 = plot(tt,ii,'Linewidth',2.5,'color','black'); p2.Color(4) = 0.6;

ax = gca;
ax.XTick = [t_u,37,67,98,t_c];
ax.XTickLabel = date([t_u,37,67,98,t_c]+1);
ax.XTickLabelRotation = 45;
xline(t_u,':','inizio Lockdown')
%ylim([0 3e5])

box on
%legend([p1,p2],'$I_{bar}$','I','Location','Best');
ylabel('casi attuali confermati');

if exist('regione','var') == 1
    title(char(regione));
else
    title('Italia')
end

fig_xlim = get(gca,'XLim');

set(gca,'FontSize',12.5);

ax_fig.XLim = fig_xlim;

if isfield(data,'regione')
    switch regione
        case 'Veneto'
            ax_fig.YLim = [0 2e4];
        case 'Lombardia'
            ax_fig.YLim = [0 6.5e4];
        case 'Emilia-Romagna'
            ax_fig.YLim = [0 3e4];
        otherwise
            ax_fig.YLim = [0 3.5e5];
    end
else
    ax_fig.YLim = [0 2e5];
end

p = get(gca, 'Position');
h = axes('Parent',gcf,'Position', [p(1)+.46 p(2)+.46 p(3)-.5 p(4)-.5],'box','on');

hold(h,'on')

%%% INSERISCO LA CURVA K

tt = linspace(t_u,t_c,50);
plot(days,K_disc,'o',...
        'MarkerSize',4,...
        'MarkerEdgeColor',[.5 .7 .1],...
        'MarkerFaceColor',[.8 .9 0]);

set(gca,'FontSize',8);

hold on
p2 = plot(tt,Kfun(tt'),'black','LineWidth',2.5);
p2.Color(4) = 0.7;
ylabel("$\kappa$",'FontSize',13)

ax = gca;
ax.XTick = [t_u,37,67,t_c];
ax.XTickLabel = date([t_u,37,67,t_c]+1);
ax.XTickLabelRotation = 45;

if ssave == 1
    if exist('regione','var') == 1
        exportgraphics(fig,...
        ['figure/' char(regione) '/riepilogo ' num2str(date(t_c+1)) '.pdf'],...
        'ContentType','vector',...
        'BackgroundColor','none')
    else
        exportgraphics(fig,...
        ['figure/Italia/riepilogo ' num2str(date(t_c+1)) '.pdf'],...
        'ContentType','vector',...
        'BackgroundColor','none')
    end
end

end