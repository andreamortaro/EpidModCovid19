% Modelli monte carlo?

clear all
close all
clc

ssave = 0;

% parametri modello
k       = 1000;
gamma   = 0.0495;
R_0     = 6.2352;
beta    = R_0*gamma;

% parametri
tstar   = 60;
nstep   = 1000;
tspan   = linspace(0,tstar,nstep);

N	= 1;
p = 40;
M = 30;
I0span_mis = 0.01;

SI = @(t,x) [-(beta - x(1)*x(2)/k)*x(1)*x(2);
              (beta - x(1)*x(2)/k)*x(1)*x(2) - gamma*x(2)];
Jac = @(t,x) [ -beta*x(2) + 2*x(1)*(x(2)^2)/k, -beta*x(1) + 2*(x(1)^2)*x(2)/k;
                beta*x(2) - 2*x(1)*(x(2)^2)/k,  beta*x(1) - 2*(x(1)^2)*x(2)/k - gamma]; 
options.Jacobian = Jac;

%% simulazione variando gli infetti

for m = I0span_mis      % run su tutti gli I0 che gli passo
    
    % salvo simulazioni
    I_hist = zeros(length(tspan),M);
    tt = zeros(length(tspan),M);
    
    fig = figure();
    hold on
    
    % imposto latex come inteprete per i grafici
    set(groot,...
        'defaulttextinterpreter','latex',...
        'defaultAxesTickLabelInterpreter','latex',...
        'defaultLegendInterpreter','latex');
        
    for j = 1:M     % fissato I0
    
        eps = rand;         % Uniformly distributed random numbers
        I0 = m*(1+p*eps);
        %I0 = m + p*eps;
        S0 = N - I0; x0 = [S0;I0]/N;
        sprintf('Iterazione: %i.\nI0 = %f, eps=%f',j,I0,eps)
        
        [t, x]  = eulerorosenbrock(SI,tspan,x0,options);
        x = x.*N;

        plot(t,x(:,2))
        hold on
        
        % salvo la simulazione
        I_hist(:,j) = x(:,2);   % lo salvo in colonna
        tt(:,j) = t;
    end
    
    grid on
    box on
    
    txt = sprintf('I0=%0.2f, M =%i, p=%.2f',m,M,p);
    title(txt)

    set(gca,'FontSize',12.5);

    if ssave == 1
        exportgraphics(fig,[txt '.pdf'],'ContentType','vector',...
                    'BackgroundColor','none')
    end
end

% media empirica
Imedio = sum(I_hist,2)./M;
h = plot(t,Imedio,'black','LineWidth',3,'Linestyle','--');
legend(h,'$\mu$')

%% Calcolo Varianza
% per ogni tempo tk calcolo la varianza, cioe' distanza media dei valori I(tk)
% dal valore medio Imedio

% welford
Var = zeros(length(tspan),1);
for tk = 1:1:length(t)
    Var(tk) = online_variance(I_hist(tk,:));
end

% confronto tra calcolo e funzione built-in
fig2 = figure();
plot(t,Var,t,var(I_hist'))
xlabel('t')
title('Varianza')
legend('welford','Var matlab')
grid on
set(gca,'FontSize',12.5);

%% plot media e varianza

figure();
hold on;      % plot media e varianza

H1 = plot(t, Imedio, 'Color', 'k', 'LineWidth', 2);
hold on
H2 = plot(t, Imedio - 2*Var,...
          t, Imedio + 2*Var,'Color', 'b');
H3 = plot(t, Imedio - 4*Var,...
          t, Imedio + 4*Var,'Color', 'b');
% H3 = plot(t, Imedio - Var,...
%           t, Imedio + Var, 'Color', 'm');
tt = [t', fliplr(t')]';
hold on
yy = [H2(1).YData, fliplr(H2(2).YData)]'; h2 = fill(tt,yy,'b','facealpha',0.2);
yy = [H3(1).YData, fliplr(H3(2).YData)]'; h3 = fill(tt,yy,'cyan','facealpha',0.2);
legend([H1, H2(1),H3(1)], '$\mu$', '2$\sigma$','4$\sigma$','Location', 'Best');
grid on
box on
axis tight

set(gca,'FontSize',12.5);
if ssave == 1
    exportgraphics(fig2,'varianza2.pdf','ContentType','vector',...
                'BackgroundColor','none')
end


% Algoritmo Welford per la varianza (varianza interattiva)
% scelto per evitare cancellazione catastrofica
function variance = online_variance(data)
n = 0;      % numero dati
media = 0;
M2 = 0;
for x = data
     n = n + 1;
     delta = x - media;
     media = media + delta/n;     % costruisco la media empirica ad ogni passo
     M2 = M2 + delta*(x - media);
     variance = M2/(n - 1);
end
end