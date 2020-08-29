
clear all
close all
clc

t0      = 0; tstar   = 20;

gamma   = 4;
R_0     = 3;
beta    = R_0*gamma;
k       = 1e-2;

N	= 1e5;

I0 = 1; S0 = N-I0;
x0 = [S0;I0]/N;          % dato iniziale in percentuale
tspan = linspace(t0,tstar,1500);


SI = @(t,x) [-(beta - x(1)*x(2)/k)*x(1)*x(2);
              (beta - x(1)*x(2)/k)*x(1)*x(2) - gamma*x(2)];
Jac = @(t,x) [ -beta*x(2) + 2*x(1)*(x(2)^2)/k, -beta*x(1) + 2*(x(1)^2)*x(2)/k;
                beta*x(2) - 2*x(1)*(x(2)^2)/k,  beta*x(1) - 2*(x(1)^2)*x(2)/k - gamma]; 
options.Jacobian = Jac;

[t, x]  = eulerorosenbrock(SI,tspan,x0,options);
x = x.*N;

epss = 0.9;
index = find([x(:,2) - epss < eps], 1, 'first');


%% figura

fig = figure();
plot(t,x(:,2),'SeriesIndex',1)

if exist('index','var')
    xline(t(index),'--','Eradication time')
    set(gca, 'XTick', unique([t(index), get(gca, 'XTick')]));
    xtickangle(-45)
end

title(['\epsilon = ' num2str(epss) ', T = ' num2str(t(index))])
xlabel("t")
axis tight
grid on
