
clear all
close all
clc

pplot = 0;      % ritardare beta

t0      = 0;
tstar   = 50;
t1      = 0;
k       = 1e3;

N	= 1e5;
I0 = 1; S0 = N-I0;

x0 = [S0;I0]/N;          % dato iniziale in percentuale
tspan = linspace(t0,tstar,1500);

gamma = 5;

%% 1.
R_0 = 2;
beta = R_0*gamma;

SI = @(t,x) [-(beta - x(1)*x(2)/k)*x(1)*x(2);
              (beta - x(1)*x(2)/k)*x(1)*x(2) - gamma*x(2)];
Jac = @(t,x) [ -beta*x(2) + 2*x(1)*(x(2)^2)/k, -beta*x(1) + 2*(x(1)^2)*x(2)/k;
                beta*x(2) - 2*x(1)*(x(2)^2)/k,  beta*x(1) - 2*(x(1)^2)*x(2)/k - gamma]; 
options.Jacobian = Jac;

[t, x]  = eulerorosenbrock(SI,tspan,x0,options);
x = x.*N;


fig = figure();
plot(t,x(:,2),'black')
hold on
area(t,x(:,2),'FaceColor','b','FaceAlpha',.3,'EdgeAlpha',.3)

%% 2.

R_0 = 1.5;
beta = R_0*gamma;


SI = @(t,x) [-(beta - x(1)*x(2)/k)*(t>=t1)*x(1)*x(2);
              (beta - x(1)*x(2)/k)*(t>=t1)*x(1)*x(2) - gamma*x(2)];

Jac = @(t,x) [ (-beta*x(2) + 2*x(1)*(x(2)^2)/k)*(t>=t1), (-beta*x(1) + 2*(x(1)^2)*x(2)/k)*(t>=t1);
                (beta*x(2) - 2*x(1)*(x(2)^2)/k)*(t>=t1), (beta*x(1) - 2*(x(1)^2)*x(2)/k)*(t>=t1) - gamma];
options.Jacobian = Jac;

tspan = linspace(t0,tstar,1000);
[t, x]  = eulerorosenbrock(SI,tspan,x0,options);
x = x.*N;

hold on
plot(t,x(:,2),'black','LineStyle','--')
area(t,x(:,2),'FaceColor','r','FaceAlpha',.3,'EdgeAlpha',.3)
%title(['beta = ' num2str(beta) ', gamma = ' num2str(gamma)])

% imposto latex come inteprete per i grafici
set(groot,...
    'defaulttextinterpreter','latex',...
    'defaultAxesTickLabelInterpreter','latex',...
    'defaultLegendInterpreter','latex');

xlabel("Tempo dal primo caso")
ylabel("Numero di casi")
text()
grid on
box off
set(gca,'XTickLabel',{},'YTickLabel',{})
% set(gca,'YTick',[], 'XTick',[])
axis([0 7 0 18000])
set(gca,'FontSize',12.5);

exportgraphics(fig,'flatteningthecurve.pdf','ContentType','vector',...
                   'BackgroundColor','none')


%% Richiesta Albi

if pplot == 1

    t0      = 0;
    tstar   = 100;
    t1      = 0.25;

    beta	=   8;   
    gamma   =   7.75;

    N	= 1e4;
    I0 = 1; S0 = N-I0;

    x0 = [S0;I0]/N;          % dato iniziale in percentuale
    tspan = linspace(t0,tstar,1000);

    %% beta sempre attivo

    SI = @(t,x) [-(beta - x(1)*x(2)/k)*x(1)*x(2);
                  (beta - x(1)*x(2)/k)*x(1)*x(2) - gamma*x(2)];

    Jac = @(t,x) [ -beta*x(2) + 2*x(1)*(x(2)^2)/k, -beta*x(1) + 2*(x(1)^2)*x(2)/k;
                    beta*x(2) - 2*x(1)*(x(2)^2)/k,  beta*x(1) - 2*(x(1)^2)*x(2)/k - gamma];
    options.Jacobian = Jac;

    [t, x]  = eulerorosenbrock(SI,tspan,x0,options);
    x = x.*N;

    fig2 = figure();
    plot(t,x(:,2),'black')
    hold on
    area(t,x(:,2),'FaceColor','b','FaceAlpha',.3,'EdgeAlpha',.3)


    %% beta attivo dopo t1

    SI = @(t,x) [-(beta - x(1)*x(2)/k)*(t>=t1)*x(1)*x(2);
                  (beta - x(1)*x(2)/k)*(t>=t1)*x(1)*x(2) - gamma*x(2)];

    Jac = @(t,x) [ (-beta*x(2) + 2*x(1)*(x(2)^2)/k)*(t>=t1), (-beta*x(1) + 2*(x(1)^2)*x(2)/k)*(t>=t1);
                    (beta*x(2) - 2*x(1)*(x(2)^2)/k)*(t>=t1), (beta*x(1) - 2*(x(1)^2)*x(2)/k)*(t>=t1) - gamma];
    options.Jacobian = Jac;


    tspan = linspace(t0,tstar,1000);
    [t, x]  = eulerorosenbrock(SI,tspan,x0,options);
    x = x.*N;

    hold on
    plot(t,x(:,2),'black','LineStyle','--')
    area(t,x(:,2),'FaceColor','r','FaceAlpha',.3,'EdgeAlpha',.3)
    hold on
    set(gca,'FontSize',12.5);
    %title(['beta = ' num2str(beta) ', gamma = ' num2str(gamma)])


    % imposto latex come inteprete per i grafici
    set(groot,...
        'defaulttextinterpreter','latex',...
        'defaultAxesTickLabelInterpreter','latex',...
        'defaultLegendInterpreter','latex');

    %xline(t1,'--','attivazione $\beta$','interpreter','latex')
    xlabel("Tempo dal primo caso")
    ylabel("Numero di casi")
    grid on
    box off
    set(gca,'XTickLabel',{},'YTickLabel',{})
    set(gca,'FontSize',12.5);
    axis([0 50 0 10])

    exportgraphics(fig2,'delaybeta.pdf','ContentType','vector',...
                       'BackgroundColor','none')
                   
end