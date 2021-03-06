function [t,x,beta,gamma] = regioni_prelock(data,K0,options)

% recupero i valori che servono
[Nass,Ibar,Rbar] = data.value;
[t_0,t_1,t_u,~,date] = data.time;

if isfield(data,'regione')
    regione = data(1).regione;
end


if nargin == 2
    ffig = 1;
    ssave = 1;
else
    if isfield(options,'ffig')
        ffig = options.ffig;
    end
    if isfield(options,'ssave')
        ssave = options.ssave;
    end 
end

pnt = 1;
if isfield(options,'pnt')
    pnt = options.pnt;
end

%% primo intervallo [t_0,t_1]

data(1).prelock = t_0;
data(2).prelock = t_1;

[t1,x1,beta1,gamma1] = stima_beta_gamma(data,K0,pnt);

%% secondo intervallo [t_1,t_u]

data(1).prelock = t_1+1;
data(2).prelock = t_u;

[t2,x2,beta2,gamma2] = stima_beta_gamma(data,K0,pnt);

%% simulazione

R0_1 = beta1/gamma1
R0_2 = beta2/gamma2

tmp1 = length(t_0:1:t_1);
tmp2 = length(t_1+1:1:t_u);
tmp3 = length(t_0:1:t_u);

beta = (beta1*tmp1 + beta2*tmp2)/tmp3;
gamma = (gamma1*tmp1 + gamma2*tmp2)/tmp3;

% update sys, nuovi parametri trovati
SI = @(t,x) [-beta*x(1)*x(2); beta*x(1)*x(2) - gamma*x(2)];

Jac = @(t,x) [-beta*x(2), -beta*x(2);
               beta*x(2), beta*x(1) - gamma];
opt.Jacobian = Jac;
opt.InitialStep = 0.001;
I0 = Ibar(t_0+1); R0 = Rbar(t_0+1); S0 = Nass-I0-R0;
x0 = [S0;I0]/Nass;                           % dato iniziale in percentuale

[t, x]  = rk4(SI,[t_0,t_u],x0,opt);       % simulazione modello

x(:,3) = ones(length(t),1) - x(:,1) - x(:,2);   % ricavo R per post-processing
x = Nass.*x;                                    % ri-normalizzo da percentuale a Nass

tm = t_0:1:t_u;

% FIGURA
if ffig == 1
    % imposto latex come inteprete per i grafici
    set(groot,...
        'defaulttextinterpreter','latex',...
        'defaultAxesTickLabelInterpreter','latex',...
        'defaultLegendInterpreter','latex');
    
    fig = figure();
    p1 = plot(tm,Ibar(tm+1),'o',...
        'MarkerSize',4,...
        'MarkerEdgeColor','red',...
        'MarkerFaceColor',[1 .6 .6]);        

    hold on
    p2 = plot(tm,Rbar(tm+1),'o',...
        'MarkerSize',4,...
        'MarkerEdgeColor',[.3 .4 .6],...
        'MarkerFaceColor',[.3 .6 .8]);
    
    p3 = plot(t,x(:,2),'SeriesIndex',2,'Linewidth',2.5);
    %p3 = plot(t,x(:,2),'color','black','Linewidth',2.5);

    hold on
    p4 = plot(t,x(:,3),'SeriesIndex',1,'Linewidth',2.5);
    
    p5 = plot(t1,x1(:,2),'SeriesIndex',4,'Linewidth',2.5); p5.Color(4) = 0.6;
    p6 = plot(t1,x1(:,3),'SeriesIndex',3,'Linewidth',2.5); p6.Color(4) = 0.6;
    p7 = plot(t2,x2(:,2),'SeriesIndex',4,'Linewidth',2.5); p7.Color(4) = 0.6;
    p8 = plot(t2,x2(:,3),'SeriesIndex',3,'Linewidth',2.5); p8.Color(4) = 0.6;
    
    p3.Color(4) = 0.6;
    p4.Color(4) = 0.6;

    ax = gca;
    ax.XTick = 0:7:14;
    ax.XTickLabel = date((0:7:14)+1);
    ax.XTickLabelRotation = 45;
    box on
    legend([p1,p2,p3,p4],'$I_{bar}$','$R_{bar}$','I','R','Location','NorthWest');
    ylabel('casi confermati');
    title(char(regione),'FontSize',13);
    set(gca,'FontSize',12.5)
    limsy=get(gca,'YLim');
    set(gca,'Ylim',[0 limsy(2)]);

    if ssave == 1
        if exist('regione','var') == 1
            exportgraphics(fig,'figure/' + regione + '/prelock.pdf','ContentType','vector',...
            'BackgroundColor','none')
        else
            exportgraphics(fig,'figure/italia_prelock.pdf','ContentType','vector',...
            'BackgroundColor','none')
        end
    end
end

return
