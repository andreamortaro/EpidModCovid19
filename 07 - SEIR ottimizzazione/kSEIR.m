function kSEIR(ssave)

%
%   la funzione kSEIR simula il modello kSEIR nel caso quattro differenti
%   tempi di attivazione del controllo k.
%

clc
close all

t_u = 14;
k = 0.25;

beta = 0.35;
gamma = 0.1;
mu = 0.2;

T = 34;
nstep = 7000;
tspan = linspace(t_u,T,nstep);

Nass = 100;
E0 = 10; I0 = 10; R0 = 0; S0 = Nass-E0-I0-R0;

x0=[S0;E0;I0]/Nass;          % dato iniziale
t1_hist = cell(0);

odeOpt = 1; % EE

%% 1. CONTROLLO SU [0,T]

%sovrascrivo i t1,t2
t1 = t_u;                % inizio controllo
t2 = T;                % fine controllo

t1_hist{1} = t1;

kSEIR = @(t,x) [-(beta - x(1)*x(3)*(t>=t1 & t<=t2)/k)*x(1)*x(3);
                 (beta - x(1)*x(3)*(t>=t1 & t<=t2)/k)*x(1)*x(3) - mu*x(2)
                  mu*x(2) - gamma*x(3)];

Jac = @(t,x) [ -beta*x(3) + 2*x(1)*(x(3)^2)*(t>=t1 & t<=t2)/k,0,...
               -beta*x(1) + 2*x(3)*(x(1)^2)*(t>=t1 & t<=t2)/k;...
                beta*x(3) - 2*x(1)*(x(3)^2)*(t>=t1 & t<=t2)/k, -mu,...
                beta*x(1) - 2*x(3)*(x(1)^2)*(t>=t1 & t<=t2)/k;...
                0, mu, -gamma];
options.Jacobian = Jac;

switch odeOpt
    case 1
        [tsol1,xsol1]= euleroesplicito(kSEIR,tspan,x0);
    case 2
        [tsol1,xsol1]= trapezi(kSEIR,tspan,x0,options);
    case 3
        [tsol1,xsol1]= eulerorosenbrock(kSEIR,tspan,x0,options);
end

%% 2. CONTROLLO SU [before peak,T]

%sovrascrivo i t1,t2
t1 = t_u+2.5;                % inizio controllo
t2 = T;                % fine controllo

t1_hist{2} = t1;

kSEIR = @(t,x) [-(beta - x(1)*x(3)*(t>=t1 & t<=t2)/k)*x(1)*x(3);
                 (beta - x(1)*x(3)*(t>=t1 & t<=t2)/k)*x(1)*x(3) - mu*x(2)
                  mu*x(2) - gamma*x(3)];

Jac = @(t,x) [ -beta*x(3) + 2*x(1)*(x(3)^2)*(t>=t1 & t<=t2)/k,0,...
               -beta*x(1) + 2*x(3)*(x(1)^2)*(t>=t1 & t<=t2)/k;...
                beta*x(3) - 2*x(1)*(x(3)^2)*(t>=t1 & t<=t2)/k, -mu,...
                beta*x(1) - 2*x(3)*(x(1)^2)*(t>=t1 & t<=t2)/k;...
                0, mu, -gamma];
options.Jacobian = Jac;

switch odeOpt
    case 1
        [tsol2,xsol2]= euleroesplicito(kSEIR,tspan,x0);
    case 2
        [tsol2,xsol2]= trapezi(kSEIR,tspan,x0,options);
    case 3
        [tsol2,xsol2]= eulerorosenbrock(kSEIR,tspan,x0,options);
end

%% 3. CONTROLLO SU [time_peak,T]

%sovrascrivo i t1,t2
t1 = t_u+5;                % inizio controllo
t2 = T;                % fine controllo

t1_hist{3} = t1;

kSEIR = @(t,x) [-(beta - x(1)*x(3)*(t>=t1 & t<=t2)/k)*x(1)*x(3);
                 (beta - x(1)*x(3)*(t>=t1 & t<=t2)/k)*x(1)*x(3) - mu*x(2)
                  mu*x(2) - gamma*x(3)];

Jac = @(t,x) [ -beta*x(3) + 2*x(1)*(x(3)^2)*(t>=t1 & t<=t2)/k,0,...
               -beta*x(1) + 2*x(3)*(x(1)^2)*(t>=t1 & t<=t2)/k;...
                beta*x(3) - 2*x(1)*(x(3)^2)*(t>=t1 & t<=t2)/k, -mu,...
                beta*x(1) - 2*x(3)*(x(1)^2)*(t>=t1 & t<=t2)/k;...
                0, mu, -gamma];
options.Jacobian = Jac;

switch odeOpt
    case 1
        [tsol3,xsol3]= euleroesplicito(kSEIR,tspan,x0);
    case 2
        [tsol3,xsol3]= trapezi(kSEIR,tspan,x0,options);
    case 3
        [tsol3,xsol3]= eulerorosenbrock(kSEIR,tspan,x0,options);
end


%% 4. CONTROLLO SU [after peak,T]

%sovrascrivo i t1,t2
t1 = t_u+10;               % inizio controllo
t2 = T;                % fine controllo

t1_hist{4} = t1;

kSEIR = @(t,x) [-(beta - x(1)*x(3)*(t>=t1 & t<=t2)/k)*x(1)*x(3);
                 (beta - x(1)*x(3)*(t>=t1 & t<=t2)/k)*x(1)*x(3) - mu*x(2)
                  mu*x(2) - gamma*x(3)];

Jac = @(t,x) [ -beta*x(3) + 2*x(1)*(x(3)^2)*(t>=t1 & t<=t2)/k,0,...
               -beta*x(1) + 2*x(3)*(x(1)^2)*(t>=t1 & t<=t2)/k;...
                beta*x(3) - 2*x(1)*(x(3)^2)*(t>=t1 & t<=t2)/k, -mu,...
                beta*x(1) - 2*x(3)*(x(1)^2)*(t>=t1 & t<=t2)/k;...
                0, mu, -gamma];
options.Jacobian = Jac;

switch odeOpt
    case 1
        [tsol4,xsol4]= euleroesplicito(kSEIR,tspan,x0);
    case 2
        [tsol4,xsol4]= trapezi(kSEIR,tspan,x0,options);
    case 3
        [tsol4,xsol4]= eulerorosenbrock(kSEIR,tspan,x0,options);
end

%% plot

tsol_all = {tsol1,tsol2,tsol3,tsol4};
xsol_all = {xsol1,xsol2,xsol3,xsol4};
fflag = zeros(length(tsol1),4);

for i =1:4
    
    tsol = tsol_all{i};
    xsol = xsol_all{i}; %xsol = xsol.*Nass;
    t1 = t1_hist{i};

    %issorted(xsol(:,1),'descend')
    
    fflag(:,i) = k - xsol(:,1).*xsol(:,2)/beta;
    
    % imposto latex come inteprete per i grafici
    set(groot,...
        'defaulttextinterpreter','latex',...
        'defaultAxesTickLabelInterpreter','latex',...
        'defaultLegendInterpreter','latex');
    
    fig = figure();
    p1 = plot(tsol+6,xsol(:,2)+xsol(:,3),'SeriesIndex',1,'LineWidth',1.5);
    hold on
    p2 = plot(tsol+6,xsol(:,2),'SeriesIndex',5,'LineWidth',1.5,'LineStyle','--');
    p3 = plot(tsol+6,xsol(:,3),'SeriesIndex',2,'LineWidth',1.5,'LineStyle','--');
    xline(t1+6,':','attivazione','LineWidth',1.5);
    H = get(fig,'CurrentAxes'); 
    set(H,'YTick',[0,.2,.4,.6,.8,1],'YTickLabel',{})
    if i == 1
        set(H,'XTick',[t1+6,25,30,35,40],'XTickLabel',{'$\tau^{*}$','','','','$t_{f}$'})
    elseif i == 2
        set(H,'XTick',[20,t1+6,25,30,35,40],'XTickLabel',{'','$\tau^{*}$','','','','$t_{f}$'})
    elseif i == 3
        set(H,'XTick',[20,t1+6,30,35,40],'XTickLabel',{'','$\tau^{*}$','','','$t_{f}$'})
    elseif i == 4
        set(H,'XTick',[20,25,t1+6,35,40],'XTickLabel',{'','','$\tau^{*}$','','$t_{f}$'})
    end
    legend([p1 p2 p3], 'E(t)+I(t)','E(t)','I(t)')
    axis([20 40 0 1]);
    xlabel("t");
    ylabel("Infetti I(t)");
    grid on
    
    set(gca,'FontSize',12.5)

    
    if ssave == 1
        exportgraphics(fig,['figure/kSEIR/subplot' num2str(i) '.pdf'],'ContentType','vector',...
                       'BackgroundColor','none');
    end
end

end
