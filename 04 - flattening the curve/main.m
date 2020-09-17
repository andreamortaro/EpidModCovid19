%% Confronto curve infetti ottenute per due valori di R_0

clear all
close all
clc

t0      = 0;
tstar   = 50;
t1      = 0;
k       = 1e3;

N	= 1e5;
I0 = 1; S0 = N-I0;

x0 = [S0;I0]/N;          % dato iniziale in percentuale
tspan = linspace(t0,tstar,1500);

gamma = 5;

%% 1. Prima curva infetti

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

% imposto latex come inteprete per i grafici
set(groot,...
    'defaulttextinterpreter','latex',...
    'defaultAxesTickLabelInterpreter','latex',...
    'defaultLegendInterpreter','latex');

plot(t,x(:,2),'black')
hold on
area(t,x(:,2),'FaceColor','r','FaceAlpha',.3,'EdgeAlpha',.3)
[m,im] = max(x(:,2));
text(t(im)*0.8,m*0.25,['$\mathcal{R}_0$ = ' num2str(R_0)],'FontSize',15,'color','black')

%% 2. seconda curva infetti: cambio il valore di R_0
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
area(t,x(:,2),'FaceColor','black','FaceAlpha',.3,'EdgeAlpha',.3)

axis([0 7 0 18000])
[m,im] = max(x(:,2));
yline(m,'--','Capacità sanitaria','FontSize',11.5)
hold on

xlabel("Tempo dalla prima infezione")
ylabel("Numero di infetti")
text()
grid off
box off
set(gca,'XTickLabel',{},'YTickLabel',{})
set(gca,'YTick',[], 'XTick',[])
text(t(im)*.85,m*0.3,['$\mathcal{R}_0$ = ' num2str(R_0)],'FontSize',15,'color','black')

set(gca,'FontSize',12.5);

exportgraphics(fig,'flatteningthecurve.pdf','ContentType','vector',...
                   'BackgroundColor','none')