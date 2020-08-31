function [t,x,beta,gamma] = regioni_prelock(data,K0,options)

global t_0 t_u date regione Ibar Rbar Nass x0 tm ym tspan pnt

% recupero i valori che servono
[Nass,Ibar,Rbar] = data.value;
[t_0,t_u,~,date] = data.time;
regione = data(1).regione;

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

tm  = t_0:1:t_u;
ym = [Ibar(tm+1),Rbar(tm+1)];

pnt = 1;
if isfield(options,'pnt')
    pnt = options.pnt;
end
nstep = pnt*t_u+1;
tspan = linspace(t_0,t_u,nstep);

I0 = Ibar(1); R0 = Rbar(1); S0 = Nass-I0-R0;
x0 = [S0;I0]/Nass;                           % dato iniziale in percentuale

% Minimizzazione: cerco beta e gamma

problem.options     = optimoptions('fmincon','Display','iter');
problem.solver      = 'fmincon';
problem.objective   = @minquad_prelock;
problem.x0 = K0;
problem.lb = [0,0];                 % impongo positivi i parametri cercati

K = fmincon(problem)

beta = K(1); gamma = K(2);  %  parametri stimati 

% update sys, nuovi parametri stimati
SI = @(t,x) [-beta*x(1)*x(2); beta*x(1)*x(2) - gamma*x(2)];

Jac = @(t,x) [-beta*x(2), -beta*x(2);
               beta*x(2), beta*x(1) - gamma];
opt.Jacobian = Jac;
opt.InitialStep = 0.001;

[t, x]  = rk4(SI,[t_0,t_u],x0,opt);

x(:,3) = ones(length(t),1) - x(:,1) - x(:,2);   % ricavo R per post-processing
x = Nass.*x;                                    % ri-normalizzo da percentuale a Nass

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
    exportgraphics(fig,'figure/' + regione + '/prelock.pdf','ContentType','vector',...
                   'BackgroundColor','none')
    end
end

return
