% Modello SIR

close all
clear all

ssave = 1;      % flag salva figura

% parameters:
beta = 0.3088;
gamma  = 0.0495;
A = [0.0066387; 54.712; 42.374];
Kfun = @(t) A(1)*exp(-((t-A(2))/A(3)).^2);

% model
SIR = @(t,x) [-beta*x(1).*x(2);...
               beta*x(1).*x(2)-gamma.*x(2)];
           
SIR_control = @(t,x) [-(beta - x(1)*x(2)/Kfun(t))*x(1)*x(2);
              (beta - x(1)*x(2)/Kfun(t))*x(1)*x(2) - gamma*x(2)];

Jac_control = @(t,x) [ -beta*x(2) + 2*x(1)*(x(2)^2)/Kfun(t), -beta*x(1) + 2*(x(1)^2)*x(2)/Kfun(t);
                beta*x(2) - 2*x(1)*(x(2)^2)/Kfun(t),  beta*x(1) - 2*(x(1)^2)*x(2)/Kfun(t) - gamma];
options2.Jacobian = Jac_control;


% discretization parameters:
tstar = 80;
options1.InitialStep = 0.1; % SIR

Nass = 1000;
I0span = [1,50,100,250,500,750,850,900];

%% figura SIR

fig1 = figure();

set(groot,...
    'defaulttextinterpreter','latex',...
    'defaultAxesTickLabelInterpreter','latex',...
    'defaultLegendInterpreter','latex');

a = linspace(0,1,Nass);
plot(a,flip(a),'SeriesIndex',2,'LineWidth',1.75)
xlabel('S')
ylabel('I')
axis([0 1 0 1])
axis('square')
hold on
xline(gamma/beta,':','\gamma/\beta','LineWidth',1,'FontSize',12,...
    'LabelOrientation','horizontal')

% % Calcolo il campo vettoriale per [0,1]x[0,1] e poi lo ridimensiono in [0,Nass]x[0,Nassass]
% xval  = linspace(0,1,15);
% yval  = linspace(0,1,15);
% vectorfield(SIR,xval,yval,0,Nass)         % campo vettoriale (ridimensionato)
% hold on

for I0 = I0span
    
    % risoluzione sistema per fissato I0
    S0 = Nass-I0;
    x0=[S0;I0]./Nass;

    % modello SIR
    [t,xsol] = rk4(SIR,[0,tstar],x0,options1);
    %xsol = xsol.*Nass;

    hold on
    plot(xsol(:,1),xsol(:,2),'SeriesIndex',1,'LineWidth',1.5)  % piano delle fasi
    drawnow
    
end

grid on
set(gca,'FontSize',12.5);

if ssave == 1
    exportgraphics(fig1,'figure/piano_fasi_SIR.pdf',...
    'ContentType','vector',...
    'BackgroundColor','none')
end
