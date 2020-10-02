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
hold on
%u = ones(1,Nass);
%plot3(gamma/beta.*u,0*u,a,'LineStyle','--','LineWidth',1.5,'Color',[0 0 0]+0.2)
%line(gamma/beta*ones(1,Nass),0*ones(1,Nass),a,'LineStyle',':','LineWidth',1.5,'Color',[0 0 0]+0.2)
%text(gamma/beta+0.15,0,0.9,'$\gamma/\beta$','FontSize',12)
%plot3(gamma/beta.*u,a,0*u,'LineStyle','--','LineWidth',1.5,'Color',[0 0 0]+0.2)
%line(gamma/beta*ones(1,Nass),a,0*ones(1,Nass),'LineStyle',':','LineWidth',1.5,'Color',[0 0 0]+0.2)
%text(gamma/beta+0.15,0.85,0,'$\gamma/\beta$','FontSize',12)

xlabel('S')
ylabel('E')
zlabel('I')
grid on
hold on

I0span = 10:100:610;
E0span = 0;

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
   
I0span = 250:250:500;
E0span = 250:250:500;

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

caz = 135; cal = 15;
view(caz,cal)
set(gca,'FontSize',12.5)

if ssave == 1
    exportgraphics(fig,'figure/piano_fasi_SEIR.pdf',...
    'ContentType','vector',...
    'BackgroundColor','none')
end