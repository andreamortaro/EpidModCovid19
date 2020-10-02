function [K,A] = stima_kcontinuo(data,K0_cont,ffig,ssave)

%
%   [K,A] = stima_kcontinuo(data,K0_cont,ffig,ssave)
%
%   Fitto i parametri K_disc che si trovano in data.
%
%   INPUTS:
%   data        : struttura contenente i dati utili
%   K0_cont     : guess iniziale
%   ffig        : flag sulla stampa della figura
%   ssave       : flag print della figura
%
%   OUTPUTS:
%   K           : funzione risultante dall'ottimizzazione, fitta i K_disc
%   A           : parametri stimati dall'ottimizzazione
%

global days K_disc Nass t_u t_c regione

% recupero i valori che servono
Nass = data(1).value;
[~,~,t_u,t_c,date] = data.time;
[days,K_disc] = data.fittingK;

if isfield(data,'regione')
    regione = data(1).regione;
    if regione == 'Lombardia'
        days = days(10:end);
        K_disc = K_disc(10:end);
    end
end


problem.options    = optimoptions('fmincon','Display','iter');
problem.solver     = 'fmincon';
problem.objective  = @minquad_kcontinuo;          % funzionale obiettivo minimizzare
problem.x0         = K0_cont;                     % guess iniziale
problem.nonlcon    = @(A)mycon(A);               % vincolo non lineare su k (=beta>0)

% %prove
problem.OptimalityTolerance = 1e-12;
problem.StepTolerance = 1e-12;
problem.FunctionTolerance = 1e-12;
problem.CostrainTolerance = 1e-2;
%problem.Algorithm = 'active-set';
problem.MaxFunctionEvaluations = 4000; %max per fmincon con punto interno 3000
problem.MaxIterations = 500;


A = fmincon(problem);

if isempty(regione) % nel caso italia in questa funzione regione = []
    K = @(t) A(1)*exp(-((t-A(2))/A(3)).^2);
else
    switch regione
        case "Veneto"
            K = @(t) A(1)*exp(-((t-A(2))/A(3)).^2).*(t<101)+A(1)*exp(-((101-A(2))/A(3)).^2).*(t>=101);
        case "Emilia-Romagna"
            K = @(t) A(1)*exp(-((t-A(2))/A(3)).^2).*(t<107)+A(1)*exp(-((107-A(2))/A(3)).^2).*(t>=107);
        otherwise
            K = @(t) A(1)*exp(-((t-A(2))/A(3)).^2);
    end
end

% if isempty(regione) % nel caso italia in questa funzione regione = []
%     K = @(t) A(1)*exp(-A(2)*t).*(1-exp(-A(3)*t)).^3;
% else
%     switch regione
%         case "Veneto"
%             K = @(t) A(1)*exp(-A(2)*t).*(1-exp(-A(3)*t)).^3.*(t<105)+A(1)*exp(-A(2)*105).*(1-exp(-A(3)*105)).^3.*(t>=105);
%         case "Emilia-Romagna"
%             K = @(t) A(1)*exp(-A(2)*t).*(1-exp(-A(3)*t)).^3.*(t<107)+A(1)*exp(-A(2)*107).*(1-exp(-A(3)*107)).^3.*(t>=107);
%         otherwise
%             K = @(t) A(1)*exp(-A(2)*t).*(1-exp(-A(3)*t)).^3;
%     end
% end


if ffig == 1
    
    % imposto latex come inteprete per i grafici
    set(groot,...
        'defaulttextinterpreter','latex',...
        'defaultAxesTickLabelInterpreter','latex',...
        'defaultLegendInterpreter','latex');

    fitting = figure();
    
    tt = linspace(t_u,t_c,50);
    plot(days,K_disc,'o',...
            'MarkerSize',4,...
            'MarkerEdgeColor',[.5 .7 .1],...
            'MarkerFaceColor',[.8 .9 0]);
    
    hold on
    p2 = plot(tt,K(tt'),'black','LineWidth',2.5);
    p2.Color(4) = 0.7;
    ylabel("$\kappa$")
    
    if exist('regione','var') == 1
        title(char(regione));
    else
        title('Italia')
    end

    ax = gca;
    ax.XTick = [t_u,37,67,98,t_c];
    ax.XTickLabel = date([t_u,37,67,98,t_c]+1);
    ax.XTickLabelRotation = 45;
    
    set(gca,'FontSize',12.5);
   
    if ssave == 1
        
        if isempty(regione) == 1
            exportgraphics(fitting,'figure/Italia/fittingk.pdf',...
            'ContentType','vector',...
            'BackgroundColor','none')
        else
            exportgraphics(fitting,'figure/' + regione + '/fittingk.pdf',...
            'ContentType','vector',...
            'BackgroundColor','none')
        end
        
    end
end

end

% VINCOLO NON LINEARE SU K

function [c,ceq] = mycon(A)

global x0 beta gamma t_u t_c regione

    %K = @(t) A(1)*exp(-A(2)*t).*(1-exp(-A(3)*t)).^3;
    K = @(t) A(1)*exp(-((t-A(2))/A(3)).^2);
    
    
    SI = @(t,x) [-(beta - x(1)*x(2)/K(t))*x(1)*x(2);
                  (beta - x(1)*x(2)/K(t))*x(1)*x(2) - gamma*x(2)];
          
    Jac = @(t,x) [ -beta*x(2) + 2*x(1)*(x(2)^2)/K(t), -beta*x(1) + 2*(x(1)^2)*x(2)/K(t);
                    beta*x(2) - 2*x(1)*(x(2)^2)/K(t),  beta*x(1) - 2*(x(1)^2)*x(2)/K(t) - gamma];
    options.Jacobian = Jac;
    
    if isempty(regione) % nel caso italia in questa funzione regione = [] 
        nstep = 100;
        tspan = linspace(t_u,90,nstep);
    else
        switch regione
            case "Veneto"
                nstep = 30;
                tspan = linspace(t_u,50,nstep);
            case "Emilia-Romagna"
                nstep = 50;
                tspan = linspace(t_u,60,nstep);
            otherwise
                nstep = 500;
                tspan = linspace(t_u,110,nstep);
        end
    end
    
    [t, xm]  = eulerorosenbrock(SI,tspan,x0,options);
    
	% Nonlinear inequality constraints (c(K)<=0)
    c = xm(:,1).*xm(:,2) - beta*K(t);   % deve essere <=0 (ATTENZIONE all'=)
    
    ceq = [];
end
