function riepilogoStat(data,tt,ii,vv,hist,options)

[t_0,~,t_u,t_c,date] = data.time;
[~,Ibar,~,totCases] = data.value;
[~,~,~,M,B] = data.parametersStat;

ssave = 1;
if isfield(options,'ssave')
    ssave = options.ssave;
end

if isfield(data,'regione')
    regione = data(1).regione;
end

[cTmedio,cTVarmedia] = casiTotaliStat(data,tt,hist);

fig = figure();
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

sm = sqrt(vv);                % scarto quadratico medio

z_alpha1 = 1.96;                    % alpha = 0.05 (confidenza) --> liv conf 0.95
delta1 = z_alpha1*sm/sqrt(M*B);     % giusto M*B?

z_alpha2 = 1.2816;                  % alpha = 0.20 (confidenza) --> liv conf 0.80
delta2 = z_alpha2*sm/sqrt(M*B);

z_alpha3 = 0.67;                    % liv confidenza 0.50
delta3 = z_alpha3*sm/sqrt(M*B);

H1 = plot(tt, ii, 'Color', [.1 .3 .5], 'LineWidth', 2);
hold on
H2 = plot(tt, ii - delta1,...
          tt, ii + delta1,'Color', 'b','Visible','Off');
H3 = plot(tt, ii - delta3,...
          tt, ii + delta3,'Color', 'b','Visible','Off');
ttt = [tt', fliplr(tt')]';
hold on
yyy = [H2(1).YData, fliplr(H2(2).YData)]'; fill(ttt,yyy,'b','facealpha',0.2);
yyy = [H3(1).YData, fliplr(H3(2).YData)]'; fill(ttt,yyy,'b','facealpha',0.5);

%legend([H1, H2(1), H3(1)], '$\mu$', '0.95\%','0.50\%','Location', 'Best');

ax = gca;
ax.XTick = [t_u,37,67,98,t_c];
ax.XTickLabel = date([t_u,37,67,98,t_c]+1);
ax.XTickLabelRotation = 45;
xline(t_u,':','inizio Lockdown')
ylabel('casi attuali attesi');

grid off
box on

if exist('regione','var') == 1
    title(char(regione));
else
    title('Italia')
end

set(gca,'FontSize',12.5);

p = get(gca, 'Position');
h = axes('Parent',gcf,'Position', [p(1)+.46 p(2)+.46 p(3)-.5 p(4)-.5],'box','on');

%% INSERISCO IL TOTALE DI CASI POSITIVI

% grafico totale casi riportati nella PC

hold(h,'on')

plot(t_0:t_c,totCases,'o',...
    'MarkerSize',4,...
    'MarkerEdgeColor',[.9 .5 .1],...
    'MarkerFaceColor',[.9 .7 .3]);
      
hold(h,'on')
plot(tt,cTmedio,'Linewidth',2.5,'color','black'); p2.Color(4) = 0.6;

% grafico con intervallo di confidenza

sm = sqrt(cTVarmedia);                % scarto quadratico medio

z_alpha1 = 1.96;                    % alpha = 0.05 (confidenza) --> liv conf 0.95
delta1 = z_alpha1*sm/sqrt(M*B);     % giusto M*B?

z_alpha2 = 1.2816;                  % alpha = 0.20 (confidenza) --> liv conf 0.80
delta2 = z_alpha2*sm/sqrt(M*B);

z_alpha3 = 0.67;                    % liv confidenza 0.50
delta3 = z_alpha3*sm/sqrt(M*B);

H1 = plot(tt, cTmedio, 'Color', [.7 0 0], 'LineWidth', 2);
hold on
H2 = plot(tt, cTmedio - delta1,...
          tt, cTmedio + delta1,'Color', 'b','Visible','Off');
H3 = plot(tt, cTmedio - delta3,...
          tt, cTmedio + delta3,'Color', 'b','Visible','Off');
ttt = [tt', fliplr(tt')]';
hold on
yyy = [H2(1).YData, fliplr(H2(2).YData)]'; fill(ttt,yyy,[.6 .8 .2],'facealpha',0.3);
yyy = [H3(1).YData, fliplr(H3(2).YData)]'; fill(ttt,yyy,[.6 .8 .2],'facealpha',0.6);

ax = gca;
ax.XTick = [t_u,37,67,98,t_c];
ax.XTickLabel = date([t_u,37,67,98,t_c]+1);
ax.XTickLabelRotation = 45;
xline(t_u,':','inizio Lockdown','FontSize',6.5)
ylabel('casi totali attesi');

grid off
box on

set(gca,'FontSize',8);

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