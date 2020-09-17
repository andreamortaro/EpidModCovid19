% imposto latex come inteprete per i grafici
set(groot,...
    'defaulttextinterpreter','latex',...
    'defaultAxesTickLabelInterpreter','latex',...
    'defaultLegendInterpreter','latex');


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
        exportgraphics(fig1,'figure/funzionale1.pdf','ContentType','vector',...
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
        exportgraphics(fig2,'figure/funzionale2.pdf','ContentType','vector',...
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
    exportgraphics(fig,'figure/funzionalezoom.pdf','ContentType','vector',...
                   'BackgroundColor','none')
end