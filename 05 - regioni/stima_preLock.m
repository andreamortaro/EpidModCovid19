%{
Stima dei parametri di interesse nel periodo precedente al Lockdown.
Cerchiamo beta e gamma tra il 24 Feb e il 09 Mar
%}

function K = stima_preLock()

global t_0 t_u date regione Ibar Rbar Nass x0 tm ym tspan pnt

tm  = t_0:1:t_u;
ym = [Ibar(tm+1),Rbar(tm+1)];

K0 = [0.5,0.05];                  % guess iniziale [beta,gamma]
pnt = 10;                         % piu nodi per migliore risoluzione sistema minquad
nstep = pnt*t_u+1;
tspan = linspace(t_0,t_u,nstep);

I0 = Ibar(1); R0 = Rbar(1); S0 = Nass-I0-R0;
x0 = [S0;I0]/Nass;                           % dato iniziale in percentuale

% Minimizzazione: cerco beta e gamma

problem.options     = optimoptions('fmincon','Display','iter');
problem.solver      = 'fmincon';
problem.objective   = @minquad_preLock;
problem.x0 = K0;
problem.lb = [0,0];                 % impongo positivi i parametri cercati

K = fmincon(problem)

k1 = K(1); k2 = K(2);  %  parametri stimati 


% update sys, nuovi parametri stimati
SI = @(t,x) [-k1*x(1)*x(2); k1*x(1)*x(2) - k2*x(2)];

Jac = @(t,x) [-k1*x(2), -k1*x(2);
               k1*x(2), k1*x(1) - k2];
options.Jacobian = Jac;

[t, x]  = rk4(SI,[t_0,t_u],x0,options);

x(:,3) = ones(length(t),1) - x(:,1) - x(:,2);   % ricavo R per post-processing
x = Nass.*x;                                    % ri-normalizzo da percentuale a Nass


% FIGURA

% imposto latex come inteprete per i grafici
set(groot,...
    'defaulttextinterpreter','latex',...
    'defaultAxesTickLabelInterpreter','latex',...
    'defaultLegendInterpreter','latex');

fig_preLock = figure();
%plot(t,x(:,2),'r',t,x(:,3),'b',tm,Ibar(tm+1),'rx',tm,Rbar(tm+1),'bx');
hold on
plot(t,x(:,2),'SeriesIndex',1);
plot(t,x(:,3),'SeriesIndex',2);
plot(tm,Ibar(tm+1),'SeriesIndex',1,'LineStyle','none','Marker','*');
plot(tm,Rbar(tm+1),'SeriesIndex',2,'LineStyle','none','Marker','*');

ax = gca;
ax.XTick = 0:7:14;
ax.XTickLabel = date((0:7:14)+1);
ax.XTickLabelRotation = 45;
%ax.YAxis.Exponent = 2;
%axis tight
box on
legend('I','R','$I_{bar}$','$R_{bar}$','Location','NorthWest');
title([char(regione), ', $\beta$ = ', num2str(k1,4), ' $\gamma$ = ', num2str(k2,3)]);
ylabel('casi confermati');
set(gca,'FontSize',12.5)


exportgraphics(fig_preLock,'figure/' + regione + '_preLock.pdf','ContentType','vector',...
               'BackgroundColor','none')

return