% Modello SIR

close all
clear all

% parameters:
alpha = 0.05; beta  = 0.5;
% model
SIR = @(t,x) [-alpha*x(1).*x(2);...
               alpha*x(1).*x(2)-beta.*x(2);...
               beta.*x(2)];

% discretization parameters:
T = 20; ntot= 100; dt = T/(ntot-1);
tspan = linspace(0,T,ntot);

S0 = 60; I0  = 40; R0 = 0;
N  = S0+I0+R0;
x0=[S0;I0;R0];
[t,xsol] = ode45(SIR,tspan,x0);

%% simulazione SIR
figure('visible','off')
plot(t,xsol(:,1),t,xsol(:,2),t,xsol(:,3))   % traiettorie di S(t),I(t),R(t)
xlabel('t')
axis('square')
title('SIR model')
legend('S','I','R')

%% diagramma fase
figure()
[Svec,Ivec] = meshgrid(1:5:N,1:5:N);
S = 1:5:N; I = 1:5:N;
vectorfield(SIR,S,I)                    % plot campo vettoriale
hold on
plot(xsol(:,1),xsol(:,2),'b','LineWidth',1.75)  % piano delle fasi
plot(1:N,N:-1:1,'LineWidth',1.75)
axis([0 N 0 N])
xlabel('S')
ylabel('I')
title('piano delle fasi')
axis('square')
drawnow
