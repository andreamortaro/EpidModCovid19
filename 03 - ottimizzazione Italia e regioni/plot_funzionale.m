% script per il plot del funzionale relativo ai k discreti

% imposto latex come inteprete per i grafici
set(groot,...
    'defaulttextinterpreter','latex',...
    'defaultAxesTickLabelInterpreter','latex',...
    'defaultLegendInterpreter','latex');

% plot con zoom
kspan = linspace(1e-12,0.05,150);
pnt = 10;
[K1,L1] = plot_minquad_kdiscreti(data,kspan,pnt);

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
[K2,L2] = plot_minquad_kdiscreti(data,kspan,pnt);
plot(K2,L2,'LineWidth',1.5)

% zoomPlot to highlight a portion of the major plot
[~,z] = zoomPlot(K1,L1,[1e-12 1e-2],[0.4 0.255 0.4 0.4],[1 3]);
z.XLim = [0 0.01];
z.YLim = [0 1e5];
z.XGrid = 'on'; z.YGrid = 'on';
z.XTick = 0:0.0025:0.010;
z.XAxis.TickLabelFormat = '%,.3f';
set(z,'FontSize',10);
hold on
legend hide

if ssave == 1
    exportgraphics(fig,'figure/Italia/funzionalezoom.pdf','ContentType','vector',...
                   'BackgroundColor','none')
end
