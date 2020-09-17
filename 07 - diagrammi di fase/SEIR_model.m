% Modello SIR

close all
clear all

ssave = 1;      % flag salva figura

% parameters:
gamma = 2;
R0 = 2.5;
beta = R0*gamma;
mu = 3;

% model
SEI = @(t,x) [-beta*x(1).*x(3);...
               beta*x(1).*x(3)-mu*x(2);...
               mu*x(2)-gamma.*x(3)];

% discretization parameters:
tstar = 5000;
options.InitialStep = 0.1; % SIR

Nass = 1000;
I0span = 10:100:610;
E0span = 0;

%% figura SEI

fig = figure();

set(groot,...
    'defaulttextinterpreter','latex',...
    'defaultAxesTickLabelInterpreter','latex',...
    'defaultLegendInterpreter','latex');

a = linspace(0,1,Nass);
b = flip(a);
z = zeros(1,Nass);

plot3(a,b,z,'SeriesIndex',2,'LineWidth',1.75);
hold on
plot3(z,a,b,'SeriesIndex',2,'LineWidth',1.75);
plot3(a,z,b,'SeriesIndex',2,'LineWidth',1.75);
plot3(a,z,z,'SeriesIndex',2,'LineWidth',1.25);
plot3(z,a,z,'SeriesIndex',2,'LineWidth',1.25);
plot3(z,z,a,'SeriesIndex',2,'LineWidth',1.25);


xlabel('S')
ylabel('E')
zlabel('I')
grid on
hold on

% % Calcolo il campo vettoriale per [0,1]x[0,1] e poi lo ridimensiono in [0,Nass]x[0,Nassass]
% xval  = linspace(0,1,15);
% yval  = linspace(0,1,15);
% vectorfield(SIR,xval,yval,0,Nass)         % campo vettoriale (ridimensionato)
% hold on

for I0 = I0span
    for E0 = E0span

        % risoluzione sistema per fissato I0
        S0 = Nass-I0-E0;
        x0=[S0;E0;I0]./Nass;

        [t,xsol] = rk4(SEI,[0,tstar],x0,options);
        %xsol = xsol.*Nass;

        % PLOT 3D
        hold on
        plot3(xsol(:,1),xsol(:,2),xsol(:,3),'SeriesIndex',1,'LineWidth',1.5)    % piano delle fasi
        xlim([0 1]); ylim([0 1]); zlim([0 1]);
        drawnow

    end
end

% X = [1; 0; 0];
% Y = [0; 1; 0];
% Z = [0; 0; 1];
 
% h2 = fill3(X,Y,Z,'r');
% h2.FaceColor=[.9 .7 .4]
% h2.FaceAlpha=0.1
% colormap(parula(5))


caz = 135; cal = 15;
%caz = 90; cal = 0;
view(caz,cal)
set(gca,'FontSize',12.5)

if ssave == 1
    exportgraphics(fig,'figure/piano_fasi_SEIR.pdf',...
    'ContentType','vector',...
    'BackgroundColor','none')
end