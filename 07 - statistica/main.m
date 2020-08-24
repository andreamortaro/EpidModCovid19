% Modelli monte carlo?

clear all
close all
clc

ssave = 0;

k       = 0.1;
gamma   = 5;
R_0     = 3.5;
beta = R_0*gamma;

tstar   = 5;
nstep = 1000;
tspan = linspace(0,tstar,nstep);
tend = 1;   % tempo finale nel plot

N	= 1;
p = 0.1;
M = 10;
I0span_mis = 0.01;

SI = @(t,x) [-(beta - x(1)*x(2)/k)*x(1)*x(2);
              (beta - x(1)*x(2)/k)*x(1)*x(2) - gamma*x(2)];
Jac = @(t,x) [ -beta*x(2) + 2*x(1)*(x(2)^2)/k, -beta*x(1) + 2*(x(1)^2)*x(2)/k;
                beta*x(2) - 2*x(1)*(x(2)^2)/k,  beta*x(1) - 2*(x(1)^2)*x(2)/k - gamma]; 
options.Jacobian = Jac;


for m = I0span_mis
    
    I_hist = zeros(length(tspan),M);
    tt = zeros(length(tspan),M);
    txt = sprintf('I0=%0.2f, M =%i, p=%.2f',m,M,p);
    fig = figure();
    hold on
    
    for j = 1:M
    
        eps = N*rand;         % Uniformly distributed random numbers
        I0 = m + p*eps;
        S0 = N - I0;
        x0 = [S0;I0]/N;
        table(I0,eps)
        
        [t, x]  = eulerorosenbrock(SI,tspan,x0,options);
        x = x.*N;

        plot(t,x(:,2))
        hold on
        I_hist(:,j) = x(:,2);   % lo salvo in colonna
        tt(:,j) = t;
    end
    
    % salvo la figura
    title(txt)
    
	% imposto latex come inteprete per i grafici
    set(groot,...
        'defaulttextinterpreter','latex',...
        'defaultAxesTickLabelInterpreter','latex',...
        'defaultLegendInterpreter','latex');
    
    grid on
    box on
    set(gca,'FontSize',12.5);
    axis([0 tend 0 N]);
    
    if ssave == 1
        exportgraphics(fig,[txt '.pdf'],'ContentType','vector',...
                    'BackgroundColor','none')
    end
end

Imedio = sum(I_hist,2)./M;
h = plot(t,Imedio,'black','LineWidth',3,'Linestyle','--');
legend(h,'$\mu$')

% calcolo varianza
% per ogni tempo tk calcolo la varianza, distanza media dei valori I(tk)
% dal valore medio Imedio
Var = zeros(length(tspan),1);
for tk = 1:1:length(t);
    tmp=0;
    for j=1:M
        tmp = tmp+ (I_hist(tk,j) - Imedio(tk))^2;
    end
    Var(tk) = tmp/M;
end

fig2 = figure();
plot(t,Var,t,var(I_hist'))
xlabel('t')
title('Varianza')
grid on
axis([0 tend 0 1e-2])
set(gca,'FontSize',12.5);
if ssave == 1
    exportgraphics(fig2,'varianza.pdf','ContentType','vector',...
                'BackgroundColor','none')
end



figure(); hold on;

mean_y = Imedio;
std_y = Var;

H1 = plot(t, Imedio, 'Color', 'k', 'LineWidth', 2);
H2 = plot(t, Imedio - 2*Var,...
          t, Imedio + 2*Var,'Color', 'b');
% H3 = plot(t, Imedio - Var,...
%           t, Imedio + Var, 'Color', 'm');
hold on
tt = [t', fliplr(t')]';
hold on
yy = [H2(1).YData, fliplr(H2(2).YData)]';
h2a = fill(tt,yy,'b');
%h2b = fill(fliplr(t),fliplr(H2(2).YData'),'b');
set(h2a,'facealpha',0.6);

legend([H1, H2(1)], ...
        '$\mu$', '2$\sigma$', ...
        'Location', 'Northwest');

grid on
box on
axis([0 tend 0 0.5]);

set(gca,'FontSize',12.5);
if ssave == 1
    exportgraphics(fig2,'varianza.pdf','ContentType','vector',...
                'BackgroundColor','none')
end