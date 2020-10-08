%% Funziona su Octave. Non provato su Matlab.
%
%   analisi sull'ammissibilità del parametro k.
%

clc
clear all
close all

ssave = 1;

% parametri
b0      = 25;
gamma   = 10;
k       = 1e-3;

Nass = 100;                 % size of population costante
I0 = 1; S0 = Nass-I0;
x0 = [S0;I0]/Nass;          % dato iniziale in percentuale

tstar = 20;                 % tempo finale

SI = @(t,x) [  -(b0 - x(1)*x(2)/k)*x(1)*x(2);
                (b0 - x(1)*x(2)/k)*x(1)*x(2) - gamma*x(2)];

Jac = @(t,x) [ -b0*x(2) + 2*x(1)*(x(2)^2)/k, -b0*x(1) + 2*(x(1)^2)*x(2)/k;
                b0*x(2) - 2*x(1)*(x(2)^2)/k,  b0*x(1) - 2*(x(1)^2)*x(2)/k - gamma];
options.Jacobian = Jac;

ntot    = 201;                          % numero di nodi
tspan   = linspace(0,tstar,ntot);
[tsol, xsol] = eulerorosenbrock (SI, tspan, x0, options);

flag = (k >= xsol(:,1).*xsol(:,2)/b0);           % flag puntuale sulla condizione sul k

xsol(:,3) = ones(length(tsol),1) - xsol(:,1) - xsol(:,2);   % ricavo R per post-processing
xsol = Nass.*xsol;                                          % ri-normalizzo da percentuale a Nass


% Risoluzione SIR

figure()
subplot(1,2,1)
plot(tsol,xsol(:,1),tsol,xsol(:,2),tsol,xsol(:,3))
legend('S','I','R')
axis([0 tstar 0 Nass])
axis("square")
xlabel('time t')
grid on
title('SIR')


subplot(1,2,2)      % campo vettoriale

xspan = linspace(0,Nass,40);
yspan = linspace(0,Nass,40);
plot(xspan,0*yspan,'k--','LineWidth',2)                 % I = 0
hold on
plot(0*xspan,yspan,'k--','LineWidth',2)                 % S = 0

% Campo vettoriale per [0,1]x[0,1] e poi ridimensiono in [0,Nass]x[0,Nass]
xval  = linspace(0,1,15);
yval  = linspace(0,1,15);
vectorfield(SI,xval,yval,0,Nass)         % cv (ridimensionato)
hold on

plot(xsol(:,1),xsol(:,2),'b','LineWidth',1.75)          % I(S)
plot(1:Nass,Nass:-1:1,'LineWidth',1.75)                 % limite dati iniziali ammissibili
plot(Nass,0,'ro','LineWidth',2)                         % equilibrio 

axis([0 Nass 0 Nass])
xlabel('S')
ylabel('I')
title('piano delle fasi')
grid on
drawnow


%% Valori ammissibili k

srange  = 0:0.01:1;
counter = 0;
krange = zeros(1,length(srange));
for s = srange
    counter=counter+1;
    krange(counter) = s*(1-s)/b0;
end


% Stampo la figura

% graphics_toolkit gnuplot
fig = figure();

set(groot,...
    'defaulttextinterpreter','latex',...
    'defaultAxesTickLabelInterpreter','latex',...
    'defaultLegendInterpreter','latex');

max_k = max(krange);
H = area(srange,[krange',2*max_k-krange'],'facecolor',[0.9 0.7 0.3]);
hold on
H2 = area(srange,krange','facecolor',"white");
%title(["Valori ammissibili per k \n\
%        fissati $\\beta_0$ = " num2str(b0) ", $\\gamma$ = " num2str(gamma)])
xlabel("$S_0$")
ylabel("$\kappa$")
axis([0 1 0 0.012])

set(gca, 'fontsize', 12.5)

if ssave == 1
    exportgraphics(fig,'figure/ammissibili.pdf',...
    'ContentType','vector',...
    'BackgroundColor','none')
end
