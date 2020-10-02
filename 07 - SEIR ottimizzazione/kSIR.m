function []=kSIR(data,ssave)

clc
close all

[Nass,Ibar,Rbar] = data.value;
[~,~,t_u,t_c,~] = data.time;
[beta,gamma] = data.parameters;
Kfun = data(3).Kvalue;

T = t_c;
tspan = linspace(t_u,t_c,50000);

Nass = 100;
% I0  = Ibar(t_u+1); R0 = Rbar(t_u+1); S0 = Nass-I0-R0;
I0 = 20; R0 = 0; S0 = Nass-I0-R0;

x0=[S0;I0]/Nass;          % dato iniziale
t1_hist = cell(0);

%% 1. CONTROLLO SU [0,T]

%sovrascrivo i t1,t2
t1 = t_u;                % inizio controllo
t2 = T;                % fine controllo

t1_hist{1} = t1;

kSIR = @(t,x) [-(beta - x(1)*x(2)*(t>=t1 & t<=t2)/Kfun(t))*x(1)*x(2);
               (beta - x(1)*x(2)*(t>=t1 & t<=t2)/Kfun(t))*x(1)*x(2) - gamma*x(2)];

Jac = @(t,x) [ -beta*x(2) + 2*x(1)*(x(2)^2)*(t>=t1 & t<=t2)/Kfun(t),...
               -beta*x(1) + 2*(x(1)^2)*x(2)*(t>=t1 & t<=t2)/Kfun(t);...
                beta*x(2) - 2*x(1)*(x(2)^2)*(t>=t1 & t<=t2)/Kfun(t),...
                beta*x(1) - 2*(x(1)^2)*x(2)*(t>=t1 & t<=t2)/Kfun(t) - gamma];
options.Jacobian = Jac;

[tsol1,xsol1]= trapezi(kSIR,tspan,x0,options);
%[tsol1,xsol1]= euleroesplicito(kSIR,tspan,x0);

%% 2. CONTROLLO SU [before peak,T]

%sovrascrivo i t1,t2
t1 = 19;                % inizio controllo
t2 = T;                % fine controllo

t1_hist{2} = t1;

kSIR = @(t,x) [-(beta - x(1)*x(2)*(t>=t1 & t<=t2)/Kfun(t))*x(1)*x(2);
               (beta - x(1)*x(2)*(t>=t1 & t<=t2)/Kfun(t))*x(1)*x(2) - gamma*x(2)];

Jac = @(t,x) [ -beta*x(2) + 2*x(1)*(x(2)^2)*(t>=t1 & t<=t2)/Kfun(t),...
               -beta*x(1) + 2*(x(1)^2)*x(2)*(t>=t1 & t<=t2)/Kfun(t);...
                beta*x(2) - 2*x(1)*(x(2)^2)*(t>=t1 & t<=t2)/Kfun(t),...
                beta*x(1) - 2*(x(1)^2)*x(2)*(t>=t1 & t<=t2)/Kfun(t) - gamma];
options.Jacobian = Jac;

[tsol2,xsol2]= euleroesplicito(kSIR,tspan,x0);

%% 3. CONTROLLO SU [time_peak,T]

%sovrascrivo i t1,t2
t1 = 24;                % inizio controllo
t2 = T;                % fine controllo

t1_hist{3} = t1;

kSIR = @(t,x) [-(beta - x(1)*x(2)*(t>=t1 & t<=t2)/Kfun(t))*x(1)*x(2);
               (beta - x(1)*x(2)*(t>=t1 & t<=t2)/Kfun(t))*x(1)*x(2) - gamma*x(2)];

Jac = @(t,x) [ -beta*x(2) + 2*x(1)*(x(2)^2)*(t>=t1 & t<=t2)/Kfun(t),...
               -beta*x(1) + 2*(x(1)^2)*x(2)*(t>=t1 & t<=t2)/Kfun(t);...
                beta*x(2) - 2*x(1)*(x(2)^2)*(t>=t1 & t<=t2)/Kfun(t),...
                beta*x(1) - 2*(x(1)^2)*x(2)*(t>=t1 & t<=t2)/Kfun(t) - gamma];
options.Jacobian = Jac;

[tsol3,xsol3]= euleroesplicito(kSIR,tspan,x0);

%% 4. CONTROLLO SU [after peak,T]

%sovrascrivo i t1,t2
t1 = 34;               % inizio controllo
t2 = T;                % fine controllo

t1_hist{4} = t1;

kSIR = @(t,x) [-(beta - x(1)*x(2)*(t>=t1 & t<=t2)/Kfun(t))*x(1)*x(2);
               (beta - x(1)*x(2)*(t>=t1 & t<=t2)/Kfun(t))*x(1)*x(2) - gamma*x(2)];

Jac = @(t,x) [ -beta*x(2) + 2*x(1)*(x(2)^2)*(t>=t1 & t<=t2)/Kfun(t),...
               -beta*x(1) + 2*(x(1)^2)*x(2)*(t>=t1 & t<=t2)/Kfun(t);...
                beta*x(2) - 2*x(1)*(x(2)^2)*(t>=t1 & t<=t2)/Kfun(t),...
                beta*x(1) - 2*(x(1)^2)*x(2)*(t>=t1 & t<=t2)/Kfun(t) - gamma];
options.Jacobian = Jac;

[tsol4,xsol4]= euleroesplicito(kSIR,tspan,x0);

%% plot

tsol_all = {tsol1,tsol2,tsol3,tsol4};
xsol_all = {xsol1,xsol2,xsol3,xsol4};

for i =1:4
    
    tsol = tsol_all{i};
    xsol = xsol_all{i}; %xsol = xsol.*Nass;
    t1 = t1_hist{i};

    issorted(xsol(:,1),'descend')
    
    % imposto latex come inteprete per i grafici
    set(groot,...
        'defaulttextinterpreter','latex',...
        'defaultAxesTickLabelInterpreter','latex',...
        'defaultLegendInterpreter','latex');
    

    fig = figure();
    %hax = plotyy(tsol, xsol(:,2),xs,ys);
    %set(hax,'xcolor','k','ycolor','k','ylim',[0,N])
    %set(hax(1), 'YTick',0:20:100)
    %ffont = 14;
    %set(hax(2), 'YTick',[0,u_max], 'YTickLabel',{'0','$u^{max}$'},'FontSize',ffont)
    plot(tsol+6,xsol(:,2),'SeriesIndex',1,'LineWidth',1.5)
    xline(t1+6,':','attivazione','LineWidth',1.5);
    H = get(fig,'CurrentAxes'); set(H,'YTickLabel',{})
    if t1 == 14
        set(H,'XTick',[t1+6,30,40,50,60],'XTickLabel',{'$\tau^{*}$','','','',''})
    elseif t1 > 14 && t1 <20
        set(H,'XTick',[20,t1+6,30,40,50,60],'XTickLabel',{'','$\tau^{*}$','','','','',''})
    elseif t1 == 24
        set(H,'XTick',[20,t1+6,40,50,60],'XTickLabel',{'','$\tau^{*}$','','',''})
    elseif t1 == 34
        set(H,'XTick',[20,30,t1+6,50,60],'XTickLabel',{'','','$\tau^{*}$','',''})
    end
    axis([20 60 0 1]);
    xlabel("t");
    ylabel("Infetti I(t)");
    grid on

    set(gca,'FontSize',12.5)
    
    if ssave == 1
        exportgraphics(fig,['figure/kSIR/subplot' num2str(i) '.pdf'],'ContentType','vector',...
                       'BackgroundColor','none');
    end
end

end
