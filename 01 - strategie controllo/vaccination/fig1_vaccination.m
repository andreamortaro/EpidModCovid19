%
%   fig1_vaccination simula il modello SIR con controllo mediante
%   vaccinazione nel caso quattro differenti tempi di attivazione
%   della strategia.
%

clc
close all
clear all

ssave = 1;

% discretizzazione temporale:
T       = 20;                       % tempo finale
ntot    = 100;                      % numero di nodi
dt      = T/(ntot-1);               % passo
tspan   = linspace(0,T,ntot);

% parametri
beta    = 0.008;
mu      = 0.08;

N   = 100;
I0  = 20; S0 = N-I0; R0 = 0;
alfa = 1;
u_max = 80;
x0=[S0;I0];          % dato iniziale
options.InitialStep = dt;
t1_hist = cell(0);

%% 1. CONTROLLO SU [0,T]

%sovrascrivo i t1,t2
t1 = 0;                % inizio controllo
t2 = T;                % fine controllo

t1_hist{1} = t1;

SIR = @(t,x) [-beta*x(1)*x(2)-alfa*x(1)*(t>=t1 & t<=t2);
               beta*x(1)*x(2)-mu*x(2)];

[tsol1,xsol1]= rk4(SIR,[0,T],x0,options);

% disegno dello scalino
scalino1 = zeros(ntot,1);
scalino1(tsol1>=t1 & tsol1<=t2) = u_max;
[xs1,ys1] = stairs(tsol1,scalino1);


%% 2. CONTROLLO SU [before peak,T]

%sovrascrivo i t1,t2
t1 = 2.5;                % inizio controllo
t2 = T;                % fine controllo

t1_hist{2} = t1;


SIR = @(t,x) [-beta*x(1)*x(2)-alfa*x(1)*(t>=t1 & t<=t2);
               beta*x(1)*x(2)-mu*x(2)];

[tsol2,xsol2]= rk4(SIR,[0,T],x0,options);

% disegno dello scalino
scalino2 = zeros(ntot,1);
scalino2(tsol2>=t1 & tsol2<=t2) = u_max;
[xs2,ys2] = stairs(tsol2,scalino2);


%% 3. CONTROLLO SU [time_peak,T]

%sovrascrivo i t1,t2
t1 = 5;                % inizio controllo
t2 = T;                % fine controllo

t1_hist{3} = t1;


SIR = @(t,x) [-beta*x(1)*x(2)-alfa*x(1)*(t>=t1 & t<=t2);
               beta*x(1)*x(2)-mu*x(2)];

[tsol3,xsol3]= rk4(SIR,[0,T],x0,options);

% disegno dello scalino
scalino3 = zeros(ntot,1);
scalino3(tsol3>=t1 & tsol3<=t2) = u_max;
[xs3,ys3] = stairs(tsol3,scalino3);


%% 4. CONTROLLO SU [after peak,T]

%sovrascrivo i t1,t2
t1 = 10;               % inizio controllo
t2 = T;                % fine controllo

t1_hist{4} = t1;


SIR = @(t,x) [-beta*x(1)*x(2)-alfa*x(1)*(t>=t1 & t<=t2);
               beta*x(1)*x(2)-mu*x(2)];

[tsol4,xsol4]= rk4(SIR,[0,T],x0,options);

% disegno dello scalino
scalino4 = zeros(ntot,1);
scalino4(tsol4>=t1 & tsol4<=t2) = u_max;
[xs4,ys4] = stairs(tsol4,scalino4);


%% plot

tsol_all = {tsol1,tsol2,tsol3,tsol4};
xsol_all = {xsol1,xsol2,xsol3,xsol4};
xs_all = {xs1,xs2,xs3,xs4};
ys_all = {ys1,ys2,ys3,ys4};

for i =1:4
    
    tsol = tsol_all{i};
    xsol = xsol_all{i};
    xs = xs_all{i};
    ys = ys_all{i};
    t1 = t1_hist{i};

    % imposto latex come inteprete per i grafici
    set(groot,...
        'defaulttextinterpreter','latex',...
        'defaultAxesTickLabelInterpreter','latex',...
        'defaultLegendInterpreter','latex');
    

    fig = figure();
    %hax = plotyy(tsol, xsol(:,2),xs,ys);
    %set(hax,'xcolor','k','ycolor','k','ylim',[0,N])
    %set(hax(1), 'YTick',0:20:100)
    %ffont = 14;
    %set(hax(2), 'YTick',[0,u_max], 'YTickLabel',{'0','$u^{max}$'},'FontSize',ffont)
    plot(tsol,xsol(:,2),'SeriesIndex',1,'LineWidth',1.5)
    xline(t1,':','attivazione','LineWidth',1.5);
    H = get(fig,'CurrentAxes'); set(H,'YTickLabel',{})
    if i == 1
        set(H,'XTick',[t1,5,10,15,T],'XTickLabel',{'$\tau^{*}$','','','','$t_{f}$'})
    elseif i == 2
        set(H,'XTick',[0,t1,5,10,15,T],'XTickLabel',{'','$\tau^{*}$','','','','$t_{f}$'})
    elseif i == 3
        set(H,'XTick',[0,t1,10,15,T],'XTickLabel',{'','$\tau^{*}$','','','$t_{f}$'})
    elseif i == 4
        set(H,'XTick',[0,5,t1,15,T],'XTickLabel',{'','','$\tau^{*}$','','$t_{f}$'})
    end
    axis([0 T 0 N]);
    xlabel("t");
    ylabel("Infetti I(t)");
    grid on

    set(gca,'FontSize',12.5)
    
    if ssave == 1
        exportgraphics(fig,['subplot' num2str(i) '.pdf'],'ContentType','vector',...
                       'BackgroundColor','none');
    end
end
