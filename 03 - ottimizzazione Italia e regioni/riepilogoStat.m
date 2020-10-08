function riepilogoStat(data,tt,ii,vv,hist,options)

%
%   riepilogoStat(data,tt,ii,vv,hist,options)
%
%   Date le simulazioni degli infetti nella struttura hist, la
%   funzione riepilogoStat stampa a video i casi attuali attesi con 
%   differenti livelli di confidenza.
%

[t_0,~,t_u,t_c,date] = data.time;
[~,Ibar,~,totCases] = data.value;
[~,~,~,M,B] = data.parametersStat;

deltatc = 30;

ssave = 1;
if isfield(options,'ssave')
    ssave = options.ssave;
end

if isfield(data,'regione')
    regione = data(1).regione;
end

casiTotaliStat(data,tt,hist,ssave);

fig = figure();
ax_fig = axes;
hold on;      % plot media e varianza

% imposto latex come inteprete per i grafici
set(groot,...
    'defaulttextinterpreter','latex',...
    'defaultAxesTickLabelInterpreter','latex',...
    'defaultLegendInterpreter','latex');
    
% grafico infetti PC

plot(t_0:t_c,Ibar,'o',...
        'MarkerSize',4,...
        'MarkerEdgeColor','red',...
        'MarkerFaceColor',[1 .6 .6]);

% grafico con intervallo di confidenza

%se = sqrt(vv)/sqrt(M*B);    % errore standard          
I_old = hist(5).sim;

% calcolo l'intervallo di confidenza al 95%
p = 0.05; CIFcn = @(x,p)prctile(x,abs([0,100]-(100-p)/2));
CI95 = zeros(length(tt),2);
tmp = cell2mat(I_old(:)');
for kk = 1:length(tt)
    x = tmp(kk,:);
    CII = CIFcn(x,100-p*100);
    CI95(kk,:) = CII;
end

% calcolo l'intervallo di confidenza al 50%
p = 0.50; CIFcn = @(x,p)prctile(x,abs([0,100]-(100-p)/2));
CI50 = zeros(length(tt),2);
tmp = cell2mat(I_old(:)');
for kk = 1:length(tt)
    x = tmp(kk,:);
    CII = CIFcn(x,100-p*100);
    CI50(kk,:) = CII;
end


hold on
H2 = plot(tt, CI95(:,1),...
         tt, CI95(:,2),'Visible','Off');
H3 = plot(tt, CI50(:,1),...
          tt, CI50(:,2),'Visible','Off');
ttt = [tt', fliplr(tt')]';
yyy = [H2(1).YData, fliplr(H2(2).YData)]'; fill(ttt,yyy,[.6 .8 1],'facealpha',0.6,'LineStyle','none');
yyy = [H3(1).YData, fliplr(H3(2).YData)]'; fill(ttt,yyy,[.4 .5 .8],'facealpha',0.6,'LineStyle','none');

plot(tt, ii, 'Color', [.1 .3 .5], 'LineWidth', 2);

%legend([H1, H2(1), H3(1)], '$\mu$', '0.95\%','0.50\%','Location', 'Best');

if isfield(data,'regione')
    switch regione
        case 'Veneto'
            ax_fig.YLim = [0 1e6];
        case 'Lombardia'
            ax_fig.YLim = [0 1e6];
        case 'Emilia-Romagna'
            ax_fig.YLim = [0 1e6];
        otherwise
            ax_fig.YLim = [0 1e6];
    end
else
    ax_fig.YLim = [0 5e6];
end
ax = gca;
ax.XTick = [t_u,37,67,98,t_c];
ax.XTickLabel = date([t_u,37,67,98,t_c]+1);
ax.XTickLabelRotation = 45;
xline(t_u,':','inizio Lockdown')
ylabel('casi attuali attesi');

grid off
legend off
box on

if exist('regione','var') == 1
    title(char(regione));
else
    title('Italia')
end

set(gca,'FontSize',12.5);

%% INSERISCO IL TOTALE DI CASI POSITIVI

% grafico totale casi riportati nella PC
%
% p = get(gca, 'Position');
% h = axes('Parent',gcf,'Position', [p(1)+.46 p(2)+.46 p(3)-.5 p(4)-.5],'box','on');
%hold(h,'on')
%
% plot(t_0:t_c,totCases,'o',...
%     'MarkerSize',4,...
%     'MarkerEdgeColor',[.9 .5 .1],...
%     'MarkerFaceColor',[.9 .7 .3]);
%       
% hold(h,'on')
% plot(tt,cTmedio,'Linewidth',2.5,'color','black'); p2.Color(4) = 0.6;
% 
% % grafico con intervallo di confidenza
% 
% sm = sqrt(cTVarmedia);                % scarto quadratico medio
% 
% p = 0.05; CIFcn = @(x,p)prctile(x,abs([0,100]-(100-p)/2));
% CI95cT = zeros(length(tt),2);
% for kk = 1:length(tt)
%     x = cT(kk,:);
%     CII = CIFcn(x,100-p*100);
%     CI95cT(kk,:) = CII;
% end
% 
% p = 0.50; CIFcn = @(x,p)prctile(x,abs([0,100]-(100-p)/2));
% CI50cT = zeros(length(tt),2);
% for kk = 1:length(tt)
%     x = cT(kk,:);
%     CII = CIFcn(x,100-p*100);
%     CI50cT(kk,:) = CII;
% end
% 
% hold on
% H2 = plot(tt, CI95cT(:,1),...
%           tt, CI95cT(:,2),'Color', [.5 .8 .3],'Visible','Off');
% H3 = plot(tt, CI50cT(:,1),...
%           tt, CI50cT(:,2),'Color', [.5 .7 .2],'Visible','Off');
% ttt = [tt', fliplr(tt')]';
% yyy = [H2(1).YData, fliplr(H2(2).YData)]'; fill(ttt,yyy,[.5 .8 .3],'facealpha',0.5);
% yyy = [H3(1).YData, fliplr(H3(2).YData)]'; fill(ttt,yyy,[.5 .7 .2],'facealpha',0.5);
% 
% H1 = plot(tt, cTmedio, 'Color', [.3 .5 0], 'LineWidth', 2.5);
% 
% ax = gca;
% ax.XTick = [t_u,37,67,98,t_c];
% ax.XTickLabel = date([t_u,37,67,98,t_c]+1);
% ax.XTickLabelRotation = 45;
% xline(t_u,':','inizio Lockdown','FontSize',6.5)
% ylabel('casi totali attesi');
% 
% grid off
% box on
% 
% set(gca,'FontSize',8);

%fig.CurrentAxes.XLim = [0 157];

if ssave == 1
    if exist('regione','var') == 1
        exportgraphics(fig,...
        ['figure/' char(regione) '/riepilogoStat ' num2str(date(t_c+1)) '.pdf'],...
        'ContentType','vector',...
        'BackgroundColor','none')
    else
        exportgraphics(fig,['figure/Italia/riepilogoStat ' num2str(date(t_c+1)) '.pdf'],...
         'ContentType','vector',...
        'BackgroundColor','none')
    end
end
